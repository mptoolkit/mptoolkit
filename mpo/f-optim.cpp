// -*- C++ -*- $Id: f-optim.cpp 1597 2015-09-05 20:16:58Z ianmcc $
//
// Optimized version of the F-matrix operator_prod
//

#include "mpo/operator_component.h"

// note: we could parallelize the construction of the G and H indices over jP

typedef std::complex<double> NumberType;
typedef LinearAlgebra::Matrix<NumberType> MatrixType;

// This version uses a Matrix* directly
struct GMatrixRef
{
   MatrixType const* F;
   MatrixType const* BH;

   GMatrixRef(MatrixType const& F_, MatrixType const& BH_) : F(&F_), BH(&BH_) {}
};

bool
operator<(GMatrixRef const& x, GMatrixRef const& y)
{
   // FIXME: should use std::less here, but in practice it surely doesn't matter
   return x.F < y.F || (x.F == y.F && x.BH < y.BH);
}

class GMatrixRefList
{
   public:
      GMatrixRefList() {}

      // we can reserve a maximum size
      void reserve(int MaxSize);

      // returns the index of a G-matrix (which is inserted if it didn't exist previously)
      int Lookup(MatrixType const& F, MatrixType const& BH);

      // do the actual evaluations (in parallel)
      void Evaluate();

      // after the evaluation, return the n'th matrix
      MatrixType const& operator[](int n) const
      {
	 return Data[n].second;
      }

   private:
      typedef std::map<GMatrixRef, int> IndexType;
      IndexType Indices;
      std::vector<std::pair<GMatrixRef, MatrixType> > Data;
};

void
GMatrixRefList::reserve(int MaxSize)
{
   Data.reserve(MaxSize);
}

int
GMatrixRefList::Lookup(MatrixType const& F, MatrixType const& BH)
{
   GMatrixRef Index(F, BH);
   IndexType::const_iterator I = Indices.find(Index);
   if (I == Indices.end())
   {
      int r = Data.size();
      Indices[Index] = r;
      Data.push_back(std::make_pair(Index, MatrixType()));
      return r;
   }
   return I->second;
}   

void
GMatrixRefList::Evaluate()
{
   for (unsigned n = 0; n < Data.size(); ++n)
   {
      Data[n].second = (*Data[n].first.F) * herm(*Data[n].first.BH);
   }
}

// array of coefficient * index into the GMatrices array
typedef std::vector<std::pair<NumberType, int> > HMatrixDescriptor;

MatrixType EvaluateH(HMatrixDescriptor const& H, GMatrixRefList const& G)
{
   CHECK(!H.empty())("H Matrix descriptor should not be empty!");

   MatrixType Result = H[0].first * G[H[0].second];
   for (unsigned n = 1; n < H.size(); ++n)
   {
      Result += H[n].first * G[H[n].second];
   }
   return Result;
}

struct ElementRec
{
   ElementRec(int sP_, int i_, HMatrixDescriptor const& H_)
      : sP(sP_), i(i_), H(H_) {}

   int sP, i;
   HMatrixDescriptor H;
};

struct OuterIndex
{
   OuterIndex(int iP_, int aP_, int jP_)
      : iP(iP_), aP(aP_), jP(jP_) {}

   bool empty() const { return Components.empty(); }

   // do the final evaluation.  Precondition: !empty()
   MatrixType Evaluate(StateComponent const& A, GMatrixRefList const& G) const;

   int iP, aP, jP;
   std::vector<ElementRec> Components;
};

MatrixType
OuterIndex::Evaluate(StateComponent const& A, GMatrixRefList const& G) const
{
   DEBUG_CHECK(!Components.empty());
   MatrixOperator::const_inner_iterator x = iterate_at(A[Components[0].sP].data(), iP, Components[0].i);
   DEBUG_CHECK(x);

   MatrixType Result = (*x) * EvaluateH(Components[0].H, G);

   for (unsigned n = 1; n < Components.size(); ++n)
   {
      x = iterate_at(A[Components[n].sP].data(), iP, Components[n].i);
      DEBUG_CHECK(x);
      Result += (*x) * EvaluateH(Components[n].H, G);
   }
   return Result;
}

// Result[a'](i',j') = M[s',s](a',a) A[s'](i',i) F[a](i,j) herm(B[s](j',j))
StateComponent
operator_prod(OperatorComponent const& M,
              StateComponent const& A, 
              StateComponent const& F,
              HermitianProxy<StateComponent> const& B)
{
   PRECONDITION_EQUAL(M.LocalBasis1(), A.LocalBasis());
   PRECONDITION_EQUAL(M.LocalBasis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), F.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), F.Basis1());
   DEBUG_PRECONDITION_EQUAL(F.Basis2(), B.base().Basis2());

   // firstly, assemble the required H matrix descriptors
   GMatrixRefList G;

   // this describes how to construct the final matrices
   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.Basis1(), A.Basis1(), B.base().Basis1());

   // iterate over all possible output matrix elements F'[aP](iP, jP)
   for (unsigned aP = 0; aP < M.Basis1().size(); ++aP)
   {
      for (unsigned iP = 0; iP < A.Basis1().size(); ++iP)
      {
	 for (unsigned jP = 0; jP < B.base().Basis1().size(); ++jP)
	 {
	    if (!is_transform_target(B.base().Basis1()[jP], M.Basis1()[aP], A.Basis1()[iP]))
	       continue;

	    // Now we have enough to construct the output descriptor
	    C.push_back(OuterIndex(iP, aP, jP));


	    // Iterate over A[sP](iP,i) (we don't care about the actual values, just the used elements)
	    for (unsigned sP = 0; sP < A.LocalBasis().size(); ++sP)
	    {
	       // We already know iP, so iterate over the row
	       for (MatrixOperator::MatrixType::data_type::value_type::const_iterator 
		       I = iterate(A[sP].data().vec()[iP]); I; ++I)
	       {
		  int i = I.index();

		  HMatrixDescriptor HMat;

		  // TODO: we already know AP
		  // Iterate over the components in M.  We already know aP, so we can iterate over
		  // the row of the matrix
		  for (OperatorComponent::data_type::data_type::value_type::const_iterator 
			  AA = iterate(M.data().vec()[aP]);
		       AA; ++AA)
		  {
		     int a = AA.index();
			
		     // Iterate over the irreducible components of M(aP,a)
		     for (SimpleRedOperator::const_iterator k = AA->begin(); k != AA->end(); ++k)
		     {
			// *k is an irreducible operator.  Iterate over the components of this operator.
			// We already know the index sP, so iterate over that row
			for (SimpleOperator::MatrixType::data_type::value_type::const_iterator
				S = iterate(k->data().vec()[sP]); S; ++S)
			{
			   int s = S.index();
			      
			   // The final index is j - we only need this if the
			   // element exists in both B.base()[s][jP, j] and
			   // F[a](i,j)

			   MatrixOperator::MatrixType::data_type::value_type::const_iterator
			      Bj = iterate(B.base()[s].data().vec()[jP]),
			      Fj = iterate(F[a].data().vec()[i]);

			   while (Bj && Fj)
			   {
			      if (Bj.index() < Fj.index())
			      {
				 ++Bj;
				 continue;
			      }
			      else if (Fj.index() < Bj.index())
			      {
				 ++Fj;
				 continue;
			      }
			      
			      int j = Bj.index();
			      DEBUG_CHECK_EQUAL(j, int(Fj.index()));

			      // now assemble the component
			      double Coeff = tensor_coefficient(B.base().Basis2()[j],
								B.base().LocalBasis()[s],
								B.base().Basis1()[jP],
								
								M.Basis2()[a],
								k->TransformsAs(),
								M.Basis1()[aP],
								
								A.Basis2()[i],
								A.LocalBasis()[sP],
								A.Basis1()[iP]);
			      
			      if (LinearAlgebra::norm_frob(Coeff) > 1E-14)
			      {
				 int GIndex = G.Lookup(*Fj, *Bj);
				 HMat.push_back(std::make_pair(Coeff * (*S), GIndex));
			      }

			      ++Bj;
			      ++Fj;
			   }
			}
		     }
		  }
		  if (!HMat.empty())
		     C.back().Components.push_back(ElementRec(sP, i, HMat));
	       }
	    }
	 }
      }
   }
   
   // do the evaluations
   G.Evaluate();

   for (unsigned n = 0; n < C.size(); ++n)
   {
      if (!C[n].empty())
      {
	 set_element(Result[C[n].aP].data(), C[n].iP, C[n].jP, C[n].Evaluate(A,G));
      }
   }

   Result.debug_check_structure();
   return Result;
}

StateComponent
contract_from_right(HermitianProxy<OperatorComponent> const& M,
		    StateComponent const& A, 
		    StateComponent const& F,
		    HermitianProxy<StateComponent> const& B)
{
   PRECONDITION_EQUAL(M.base().LocalBasis1(), A.LocalBasis());
   PRECONDITION_EQUAL(M.base().LocalBasis2(), B.base().LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.base().Basis2(), F.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(A.Basis2(), F.Basis1());
   DEBUG_PRECONDITION_EQUAL(F.Basis2(), B.base().Basis2());

   // firstly, assemble the required H matrix descriptors
   GMatrixRefList G;

   // this describes how to construct the final matrices
   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.base().Basis1(), A.Basis1(), B.base().Basis1());

   // iterate over all possible output matrix elements F'[aP](iP, jP)
   for (unsigned aP = 0; aP < M.base().Basis1().size(); ++aP)
   {
      for (unsigned iP = 0; iP < A.Basis1().size(); ++iP)
      {
	 for (unsigned jP = 0; jP < B.base().Basis1().size(); ++jP)
	 {
	    if (!is_transform_target(B.base().Basis1()[jP], M.base().Basis1()[aP], A.Basis1()[iP]))
	       continue;

	    // Now we have enough to construct the output descriptor
	    C.push_back(OuterIndex(iP, aP, jP));


	    // Iterate over A[sP](iP,i) (we don't care about the actual values, just the used elements)
	    for (unsigned sP = 0; sP < A.LocalBasis().size(); ++sP)
	    {
	       // We already know iP, so iterate over the row
	       for (MatrixOperator::MatrixType::data_type::value_type::const_iterator 
		       I = iterate(A[sP].data().vec()[iP]); I; ++I)
	       {
		  int i = I.index();

		  HMatrixDescriptor HMat;

		  // TODO: we already know AP
		  // Iterate over the components in M.  We already know aP, so we can iterate over
		  // the row of the matrix
		  for (OperatorComponent::data_type::data_type::value_type::const_iterator 
			  AA = iterate(M.base().data().vec()[aP]);
		       AA; ++AA)
		  {
		     int a = AA.index();
			
		     // Iterate over the irreducible components of M(aP,a)
		     for (SimpleRedOperator::const_iterator k = AA->begin(); k != AA->end(); ++k)
		     {
			// *k is an irreducible operator.  Iterate over the components of this operator.
			// We already know the index sP, so iterate over that row
			for (SimpleOperator::MatrixType::data_type::value_type::const_iterator
				S = iterate(k->data().vec()[sP]); S; ++S)
			{
			   int s = S.index();
			      
			   // The final index is j - we only need this if the
			   // element exists in both B.base()[s][jP, j] and
			   // F[a](i,j)

			   MatrixOperator::MatrixType::data_type::value_type::const_iterator
			      Bj = iterate(B.base()[s].data().vec()[jP]),
			      Fj = iterate(F[a].data().vec()[i]);

			   while (Bj && Fj)
			   {
			      if (Bj.index() < Fj.index())
			      {
				 ++Bj;
				 continue;
			      }
			      else if (Fj.index() < Bj.index())
			      {
				 ++Fj;
				 continue;
			      }
			      
			      int j = Bj.index();
			      DEBUG_CHECK_EQUAL(j, int(Fj.index()));

			      // now assemble the component
			      double Coeff = tensor_coefficient(B.base().Basis2()[j],
								B.base().LocalBasis()[s],
								B.base().Basis1()[jP],
								
								M.base().Basis2()[a],
								k->TransformsAs(),
								M.base().Basis1()[aP],
								
								A.Basis2()[i],
								A.LocalBasis()[sP],
								A.Basis1()[iP]);
			      
			      if (LinearAlgebra::norm_frob(Coeff) > 1E-14)
			      {
				 int GIndex = G.Lookup(*Fj, *Bj);
				 HMat.push_back(std::make_pair(Coeff * herm(*S), GIndex));
			      }

			      ++Bj;
			      ++Fj;
			   }
			}
		     }
		  }
		  if (!HMat.empty())
		     C.back().Components.push_back(ElementRec(sP, i, HMat));
	       }
	    }
	 }
      }
   }
   
   // do the evaluations
   G.Evaluate();

   for (unsigned n = 0; n < C.size(); ++n)
   {
      if (!C[n].empty())
      {
	 set_element(Result[C[n].aP].data(), C[n].iP, C[n].jP, C[n].Evaluate(A,G));
      }
   }

   Result.debug_check_structure();
   return Result;
}
