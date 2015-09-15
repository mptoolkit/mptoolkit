// -*- C++ -*- $Id$

//
// DensityMatrix<MatrixOperator>
//

#include "linearalgebra/matrix_utility.h"

template <typename FwdIterX>
MatrixOperator DensityMatrix<MatrixOperator>::ConstructTruncator(FwdIterX Start, FwdIterX End) const
{
   VectorBasis NewBasis(B.GetSymmetryList());

   // make a pass over the eigenvalue list and get the linear indices of the
   // states to keep for each quantum number
   int NumQ = B.size();
   std::vector<std::set<int> > LinearMapping(NumQ);
   for (FwdIterX Iter = Start; Iter != End; ++Iter)
   {
      LinearMapping[Iter->Subspace].insert(Iter->Index);
   }

   // Now we construct the truncated basis
   std::vector<int> NewSubspace(NumQ, -1); // maps the q to the new label
   for (int q = 0; q < NumQ; ++q)
   {
      if (LinearMapping[q].size() > 0)
      {
	 NewSubspace[q] = NewBasis.size();
         NewBasis.push_back(B[q], LinearMapping[q].size());
      }
   }
   // Now we construct the actual transform
   MatrixOperator Transform(NewBasis, B.MappedBasis(), QuantumNumber(B.GetSymmetryList()));
   for (std::size_t s = 0; s < Transform.Basis2().size(); ++s)
   {
      int q;
      LinearAlgebra::Range LinRange;
      boost::tie(q, LinRange) = B.Lookup(s);

      // sp is the subspace in NewBasis
      int sp = NewSubspace[q];
      if (sp == -1) continue;

      Transform(sp, s) = RawDMList[q](std::vector<int>(LinearMapping[q].begin(), LinearMapping[q].end()), LinRange);
   }
   return Transform;
}

//
// DensityMatrix<SimpleOperator>
//


template <class FwdIter>
SimpleOperator DensityMatrix<SimpleOperator>::ConstructTruncator(FwdIter Start, FwdIter End) const
{
   BasisList NewBasis(B.GetSymmetryList());
   // make a pass over the eigenvalue list and get the linear indices of the
   // states to keep for each quantum number
   std::vector<std::pair<int, int> > KeptStates;
   for (FwdIter Iter = Start; Iter != End; ++Iter)
   {
      NewBasis.push_back(B[Iter->Subspace]);
      KeptStates.push_back(std::make_pair(Iter->Subspace, Iter->Index));
   }
   DEBUG_CHECK_EQUAL(NewBasis.size(), KeptStates.size());

   // Now we construct the actual transform
   SimpleOperator Transform(NewBasis, B.MappedBasis(), QuantumNumber(B.GetSymmetryList()));
   for (std::size_t s = 0; s < Transform.Basis2().size(); ++s)
   {
      int q, qi;
      boost::tie(q, qi) = B.Lookup(s);

      for (std::size_t sp = 0; sp < NewBasis.size(); ++sp)
      {
	 int qp, qpi;
	 boost::tie(qp, qpi) = KeptStates[sp];
	 if (qp == q)
	    Transform(sp, s) = RawDMList[q](qpi, qi);
      }
   }
   return Transform;
}

template <class FwdIter>
SimpleOperator DensityMatrix<SimpleOperator>::ConstructUnnormalizedTruncator(FwdIter Start, 
									     FwdIter End) const
{
   BasisList NewBasis(B.GetSymmetryList());
   // make a pass over the eigenvalue list and get the linear indices of the
   // states to keep for each quantum number
   std::vector<std::pair<int, int> > KeptStates;
   std::vector<double> KeptEigenvalue;
   for (FwdIter Iter = Start; Iter != End; ++Iter)
   {
      NewBasis.push_back(B[Iter->Subspace]);
      KeptStates.push_back(std::make_pair(Iter->Subspace, Iter->Index));
      KeptEigenvalue.push_back(1.0 / std::sqrt(Iter->Eigenvalue));
   }
   DEBUG_CHECK_EQUAL(NewBasis.size(), KeptStates.size());

   // Now we construct the actual transform
   SimpleOperator Transform(NewBasis, B.MappedBasis(), QuantumNumber(B.GetSymmetryList()));
   for (std::size_t s = 0; s < Transform.Basis2().size(); ++s)
   {
      int q, qi;
      boost::tie(q, qi) = B.Lookup(s);

      for (std::size_t sp = 0; sp < NewBasis.size(); ++sp)
      {
	 int qp, qpi;
	 boost::tie(qp, qpi) = KeptStates[sp];
	 if (qp == q)
	    Transform(sp, s) = RawDMList[q](qpi, qi) * KeptEigenvalue[sp];
      }
   }
   return Transform;
}

//
// Singular value decomposition
//

template <typename FwdIter>
void
SingularDecomposition<StateComponent, StateComponent>::
ConstructMatrices(FwdIter first, FwdIter last, 
		  StateComponent& A, RealDiagonalOperator& C, StateComponent& B)
{
   // make a pass over the eigenvalue list and get the linear indices of the
   // states to keep for each quantum number
   int NumQ = UsedQuantumNumbers.size();
   std::vector<std::set<int> > LinearMapping(NumQ);
   for (FwdIter Iter = first; Iter != last; ++Iter)
   {
      LinearMapping[Iter->Subspace].insert(Iter->Index);
   }

   this->ConstructOrthoMatrices(LinearMapping, A, C, B);
}

template <typename FwdIter>
void
SingularDecomposition<MatrixOperator, MatrixOperator>::
ConstructMatrices(FwdIter first, FwdIter last, 
		       MatrixOperator& A, 
		       RealDiagonalOperator& C, 
		       MatrixOperator& B)
{
   // make a pass over the eigenvalue list and get the linear indices of the
   // states to keep for each quantum number
   int NumQ = UsedQuantumNumbers.size();
   std::vector<std::set<int> > LinearMapping(NumQ);
   for (FwdIter Iter = first; Iter != last; ++Iter)
   {
      LinearMapping[Iter->Subspace].insert(Iter->Index);
   }

   this->ConstructOrthoMatrices(LinearMapping, A, C, B);
}
