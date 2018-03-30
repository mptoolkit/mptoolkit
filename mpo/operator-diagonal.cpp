// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/operator-diagonal.cpp
//
// Copyright (C) 2015-2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
//

#include "mpo/operator_component.h"

// note: we could parallelize the construction of the G and H indices over jP

#if 0

// This version uses a Matrix* directly
struct JMatrixRef
{
   Matrix const* AH;
   Matrix const* E;

   JMatrixRef(Matrix const& AH_, Matrix const& E_) : AH(&AH_), E(&E_) {}
};

bool
operator<(JMatrixRef const& x, JMatrixRef const& y)
{
   // FIXME: should use std::less here, but in practice it surely doesn't matter
   return x.AH < y.AH || (x.AH == y.AH && x.E < y.E);
}

class JMatrixRefList
{
   public:
      JMatrixRefList() {}

      // we can reserve a maximum size
      void reserve(int MaxSize);

      // returns the index of a G-matrix (which is inserted if it didn't exist previously)
      int Lookup(Matrix const& AH, Matrix const& E);

      // do the actual evaluations (in parallel)
      void Evaluate();

      // after the evaluation, return the n'th matrix
      Matrix const& operator[](int n) const
      {
         return Data[n].second;
      }

   private:
      typedef std::map<JMatrixRef, int> IndexType;
      IndexType Indices;
      std::vector<std::pair<JMatrixRef, Matrix> > Data;
};

void
JMatrixRefList::reserve(int MaxSize)
{
   Data.reserve(MaxSize);
}

int
JMatrixRefList::Lookup(Matrix const& AH, Matrix const& E)
{
   JMatrixRef Index(AH, E);
   IndexType::const_iterator I = Indices.find(Index);
   if (I == Indices.end())
   {
      int r = Data.size();
      Indices[Index] = r;
      Data.push_back(std::make_pair(Index, Matrix()));
      return r;
   }
   return I->second;
}

void
JMatrixRefList::Evaluate()
{
   for (unsigned n = 0; n < Data.size(); ++n)
   {
      Data[n].second = herm(*Data[n].first.AH) * (*Data[n].first.E);
   }
}

// array of coefficient * index into the GMatrices array
typedef std::vector<std::pair<complex, int> > KMatrixDescriptor;

Matrix EvaluateK(KMatrixDescriptor const& K, JMatrixRefList const& J)
{
   CHECK(!K.empty())("K Matrix descriptor should not be empty!");

   Matrix Result = K[0].first * J[K[0].second];
   for (unsigned n = 1; n < K.size(); ++n)
   {
      Result += K[n].first * J[K[n].second];
   }
   return Result;
}

struct ElementRec
{
   ElementRec(int s_, int jP_, KMatrixDescriptor const& K_)
      : s(s_), jP(jP_), K(K_) {}

   int s, jP;
   KMatrixDescriptor K;
};

struct OuterIndex
{
   OuterIndex(int i_, int a_, int j_)
      : i(i_), a(a_), j(j_) {}

   bool empty() const { return Components.empty(); }

   // do the final evaluation.  Precondition: !empty()
   Matrix Evaluate(StateComponent const& A, JMatrixRefList const& J) const;

   int i, a, j;
   std::vector<ElementRec> Components;
};

Matrix
OuterIndex::Evaluate(StateComponent const& B, JMatrixRefList const& J) const
{
   DEBUG_CHECK(!Components.empty());
   MatrixOperator::const_inner_iterator x = iterate_at(B[Components[0].s].data(), Components[0].jP, j);
   DEBUG_CHECK(x);

   Matrix Result = EvaluateK(Components[0].K, J) * (*x);

   for (unsigned n = 1; n < Components.size(); ++n)
   {
      x = iterate_at(B[Components[n].s].data(), Components[n].jP, j);
      DEBUG_CHECK(x);
      Result += EvaluateK(Components[n].K, J) * (*x);
   }
   return Result;
}

#endif

Matrix
outer_diagonal_conj(Matrix const& M1, Matrix const& M2)
{
   using std::conj;
   Matrix Result(M1.rows(), M2.rows());
#if 0
   for (unsigned i = 0; i < M1.rows(); ++i)
   {
      for (unsigned j = 0; j < M2.cols(); ++j)
      {
	 Result(i,j) = M1(i,i) * conj(M2(j,j));
      }
   }
#endif
   return Result;
}

// Result[s](j',j) = M(s,s)[a',a] E[a'](j',j') B[s](j',j) herm(F[a](j,j))
StateComponent
operator_prod_inner_diagonal(OperatorComponent const& M,
			     StateComponent const& E,
			     HermitianProxy<StateComponent> const& F)
{
   DEBUG_PRECONDITION_EQUAL(M.Basis1(), E.LocalBasis());
   DEBUG_PRECONDITION_EQUAL(M.Basis2(), F.base().LocalBasis());

   using ::norm_frob;

   // firstly, assemble the required H matrix descriptors
   //   JMatrixRefList JList;

   // this describes how to construct the final matrices
   //   std::vector<OuterIndex> C;

   // Firstly, iterate over the outer indices - we want to be careful to only
   // add indices that are definitely needed, the MPS might have more sparseness than
   // the minimum imposed by the symmetries.
   StateComponent Result(M.LocalBasis1(), E.Basis1(), F.base().Basis1());

   // iterate over all possible output matrix elements A'[s'](i, j)
   for (unsigned sP = 0; sP < M.LocalBasis1().size(); ++sP)
   {
      for (unsigned iP = 0; iP < E.Basis1().size(); ++iP)
      {
         //double degree_i = degree(E.base().Basis1()[i]);

         for (unsigned i = 0; i < F.base().Basis1().size(); ++i)
         {
            if (!is_transform_target(F.base().Basis1()[i], M.LocalBasis1()[sP], E.Basis1()[iP]))
               continue;


	    // Now we have enough to construct the output descriptor
	    //      C.push_back(OuterIndex(iP, sP, i));

	    // Iterate over F[a](i,j)
	    for (unsigned a = 0; a < F.base().LocalBasis().size(); ++a)
	    {
	       int j = i;
	       int jP = iP;

               auto J = F.base()[a].row(i).find(j);
               if (J == F.base()[a].row(i).end())
		  continue;


	       // Iterate over the components in M[a'a]
	       // Unlike the F-matrix case, we don't know the leading index so we need to iterate
	       // over both indices
	       for (unsigned aP = 0; aP < M.Basis1().size(); ++aP)
	       {
                  auto AA = M.row(aP).find(a);
                  if (AA == M.row(aP).end())
		     continue;

		  // Iterate over the irreducible components of M(aP,a)
                  for (auto const& k : (*AA).value)
                  {
		     // k is an irreducible operator.  Iterate over the components of this operator.
		     // We already know sP
		     // We already know the index sP, so iterate over that row
                     for (auto const& S : k.row(sP))
                     {
			int s = S.col();
			if (s != sP)
			   continue;

			// The final index is j' - we only need this if the
			// element exists in both E.base()[a'](i',j')
			//and B.base()[s](j',j)

                        auto EjP = E[aP].row(iP).find(jP);
			if (EjP == E[aP].row(iP).end())
			   continue;

			//double degree_iP = degree(A.base().Basis1()[iP]);

			// now assemble the component
			double Coeff = tensor_coefficient(F.base().Basis1()[j],
							  M.LocalBasis2()[s],
							  E.Basis1()[jP],

							  M.Basis2()[a],
							  k.TransformsAs(),
							  M.Basis1()[aP],

							  F.base().Basis1()[i],
							  M.LocalBasis1()[sP],
							  E.Basis1()[iP]);

			if (norm_frob(Coeff) > 1E-14)
			{
                           Result[sP].add(iP, i, Coeff * S.value * outer_diagonal_conj(*EjP, *J));
			}
		     }
                  }
                  //              if (!KMat.empty())
                  //                 C.back().Components.push_back(ElementRec(s, jP, KMat));
               }
            }
         }
      }
   }

#if 0
   // do the evaluations
   JList.Evaluate();

   for (unsigned n = 0; n < C.size(); ++n)
   {
      if (!C[n].empty())
      {
         set_element(Result[C[n].a].data(), C[n].i, C[n].j, C[n].Evaluate(B,JList));
      }
   }
#endif

   Result.debug_check_structure();
   return Result;
}
