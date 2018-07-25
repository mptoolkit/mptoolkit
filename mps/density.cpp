// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mps/density.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "density.h"
#include "common/proccontrol.h"
#include "blas/matrix-eigen.h"


LinearBasis<VectorBasis>::LinearBasis(VectorBasis const& B)
   : VectorBasis(B.GetSymmetryList()), Orig(B), Mapping(B.size())
{
   // Calculate the total dimension of each quantum number in B, and
   // fill in the Range component of the Mapping.
   std::map<QuantumNumber, int> QNDimension;
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      int Dim = B.dim(i);
      QuantumNumber q = B[i];

      int& MappedDim(QNDimension[q]);
      Mapping[i].second = blas::Range(MappedDim, MappedDim + Dim);
      MappedDim += Dim;
   }

   // Now we know the dimensions of the quantum numbers, we can add the subspaces the basis
   std::map<QuantumNumber, int> QNSubspace;
   for (std::map<QuantumNumber, int>::const_iterator I = QNDimension.begin();
        I != QNDimension.end(); ++I)
   {
      QNSubspace[I->first] = this->size();
      this->push_back(I->first, I->second);
   }

   // finally, fill in the subspace component of Mapping
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Mapping[i].first = QNSubspace[B[i]];
   }
}

LinearBasis<BasisList>::LinearBasis(BasisList const& B)
   : VectorBasis(B.GetSymmetryList()), Orig(B), Mapping(B.size())
{
   // Calculate the total dimension of each quantum number in B, and
   // fill in the Range component of the Mapping.
   std::map<QuantumNumber, int> QNDimension;
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Mapping[i].second = QNDimension[B[i]]++;

   }

   // Now we know the dimensions of the quantum numbers, we can add the subspaces the basis
   std::map<QuantumNumber, int> QNSubspace;
   for (std::map<QuantumNumber, int>::const_iterator I = QNDimension.begin();
        I != QNDimension.end(); ++I)
   {
      QNSubspace[I->first] = this->size();
      this->push_back(I->first, I->second);
   }

   // finally, fill in the subspace component of Mapping
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Mapping[i].first = QNSubspace[B[i]];
   }
}

int
LinearBasis<BasisList>::ReverseLookup(int s, int index) const
{
   for (std::size_t i = 0; i < Mapping.size(); ++i)
   {
      if (Mapping[i].first == s && Mapping[i].second == index) return i;
   }
   return -1;
}

std::ostream&
DensityMatrixBase::DensityMatrixReport(std::ostream& outstream, int MaxEigenvalues, bool Base2, bool ShowDegen, bool Quiet)
{
   std::ostringstream out;
   out.precision(12);
   out << std::scientific;
   if (MaxEigenvalues < 0) MaxEigenvalues = EigenInfoList.size();
   if (!Quiet)
   {
      out << "#Eigenvalue sum = " << this->EigenSum() << '\n';
      out << "#von Neumann Entropy " << (Base2 ? "(base 2)" : "(base e)") << " = " << this->Entropy() << '\n';
      if (MaxEigenvalues > 0)
	 out << "#Number    #Eigenvalue         #Degen    #Weight               #Energy             #QuantumNumber\n";
   }
   int n = 0;
   int TotalDegree = 0;
   double EShift = 0;
   double EScale = 1.0;
#if 0
   if (EigenInfoList.size() >= 2)
   {
      EShift = -log(EigenInfoList[0].Eigenvalue / this->EigenSum());
      double E2 = -log(EigenInfoList[1].Eigenvalue / this->EigenSum());
      EScale = E2 - EShift;
   }
#endif
   for (const_iterator Iter = this->begin(); Iter != this->end() && n < MaxEigenvalues; ++Iter)
   {
      double EVal = Iter->Eigenvalue / this->EigenSum();
      double Energy = EVal > 0 ? ((-log(EVal) - EShift) / EScale) : 0.0;

      int OuterDegen = ShowDegen ? Iter->Degree : 1;
      int DisplayDegen = ShowDegen ? 1 : Iter->Degree;

      for (int i = 0; i < OuterDegen; ++i)
      {
         double Weight = EVal * Iter->Degree / OuterDegen;
         TotalDegree += Iter->Degree;
         ++n;
         out << std::right << std::setw(7) << n << "  "
             << std::right << std::setw(20) << Iter->Eigenvalue
             << "  " << std::setw(6) << DisplayDegen
             << "  " << std::setw(20) << Weight
             << "  " << std::setw(20) << Energy
             << "  " << std::left
             << this->Lookup(Iter->Subspace) << '\n';
      }
   }
   if (!Quiet)
   {
      out << '#' << n << " out of " << (ShowDegen ? TotalDegree : EigenInfoList.size()) << " eigenvalues shown.  ";
      out << "Total degree = " << TotalDegree << '\n';
   }
   outstream << out.str();
   return outstream;
}

double
DensityMatrixBase::Entropy(bool Base2) const
{
   double x = 0;
   for (const_iterator Iter = begin(); Iter != end(); ++Iter)
   {
      double EVal = Iter->Eigenvalue / this->EigenSum();
      if (EVal > 0)
         x -= EVal * (Base2 ? log2(EVal) : log(EVal)) * Iter->Degree;
   }
   return x;
}

double
DensityMatrixBase::EvaluateCasimir(int n) const
{
   double x = 0;
   double ESum = this->EigenSum();
   for (const_iterator Iter = begin(); Iter != end(); ++Iter)
   {
      x += (Iter->Eigenvalue / ESum) * Iter->Degree * casimir(this->Lookup(Iter->Subspace), n);
   }
   return x;
}

double
DensityMatrixBase::EvaluateCasimirMoment(int n) const
{
   double x = 0;
   double ESum = this->EigenSum();
   double c = this->EvaluateCasimir(n);
   for (const_iterator Iter = begin(); Iter != end(); ++Iter)
   {
      double xx = casimir(this->Lookup(Iter->Subspace), n);
      x += (Iter->Eigenvalue / ESum) * Iter->Degree * (xx-c) * (xx-c);
   }
   return x;
}

void DensityMatrixBase::DiagonalizeDMHelper(bool Sort)
{
   //   double Point1 = ProcControl::GetCPUTime();  // Point1 is between construction & diagonalization

   ESum = 0;  // running sum of the eigenvalues
   // diagonalize the DM
   blas::Vector<double> Eigenvalues(MaxLinearDimension);
   for (std::size_t q1 = 0; q1 < RawDMList.size(); ++q1)
   {
      int CurrentDegree = degree(this->Lookup(q1));
      //      std::cout << "Raw DM is\n" << RawDMList[q1] << std::endl;
      //      TRACE(RawDMList[q1])(RawDMList[q1].size1())(RawDMList[q1].size2());
      blas::Vector<double, blas::gpu_tag> EVal(RawDMList[q1].rows());
      DiagonalizeHermitian(RawDMList[q1], EVal);
      blas::Vector<double> Eigenvalues = get_wait(EVal);
      // add the eigenvalues and eigenvector pointers to EigenInfoList
      for (std::size_t i = 0; i < RawDMList[q1].cols(); ++i)
      {
         EigenInfoList.push_back(EigenInfo(Eigenvalues[i], CurrentDegree, q1, i));
         ESum += Eigenvalues[i] * CurrentDegree;
      }
   }

   // sort the eigenvalues
   if (Sort)
      std::sort(EigenInfoList.begin(), EigenInfoList.end(), EigenCompare);
}


//
// DensityMatrix<MatrixOperator>
//

DensityMatrix<MatrixOperator>::DensityMatrix(MatrixOperator const& Op)
   : B(Op.Basis1())
{
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), Op.Basis2());
   // initialize the RawDMList.  At the same time, get the maximum linear size of the subspaces
   RawDMList.clear();
   RawDMList.reserve(B.size());
   MaxLinearDimension = 0;
   for (std::size_t q = 0; q < B.size(); ++q)
   {
      int Dim = B.dim(q);
      //      std::cout << "dimension of " << q1 << " is " << Dim << std::endl;
      RawDMList.emplace_back(Dim, Dim, 0.0);
      MaxLinearDimension = std::max<int>(MaxLinearDimension, Dim);
   }

   // Fill the raw density matrices
   for (auto const& rOp : Op)
   {
      for (auto const& cOp : rOp)
      {
         int tp;
         blas::Range rtp;
         std::tie(tp, rtp) = B.Lookup(rOp.row());
         int t;
         blas::Range rt;
         std::tie(t, rt) = B.Lookup(cOp.col());
         CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal")(B[tp])(B[t])(Op)(B);
         RawDMList[tp](rtp, rt) = cOp.value;
      }
   }
   // diagonalize them
   this->DiagonalizeDMHelper();
}

#if 0
DensityMatrix<MatrixOperator>::DensityMatrix(MatrixOperator const& Op, MatrixOperator const& WavefunctionDM)
   : B(Op.Basis1())
{
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), Op.Basis2());
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), WavefunctionDM.Basis1());
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), WavefunctionDM.Basis2());

   // We have a second RawDMType, for the WavefunctionDM
   std::vector<RawDMType> PsiDMList(B.size());

   // initialize the RawDMList.  At the same time, get the maximum linear size of the subspaces
   RawDMList.resize(B.size());
   MaxLinearDimension = 0;
   for (std::size_t q = 0; q < B.size(); ++q)
   {
      int Dim = B.dim(q);
      //      std::cout << "dimension of " << q1 << " is " << Dim << std::endl;
      RawDMList[q] = RawDMType(Dim, Dim, 0.0);
      PsiDMList[q] = RawDMType(Dim, Dim, 0.0);
      MaxLinearDimension = std::max<int>(MaxLinearDimension, Dim);
   }

   // Fill the raw density matrices
   for (auto const& rOp : Op)
   {
      for (auto const& cOp : rOp)
      {
         int tp;
         LinearAlgebra::Range rtp;
         std::tie(tp, rtp) = B.Lookup(rOp.row());
         int t;
         LinearAlgebra::Range rt;
         std::tie(t, rt) = B.Lookup(cOp.col());
         CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal")(B[tp])(B[t])(Op)(B);
         RawDMList[tp](rtp, rt) = *J;
      }
   }

   for (auto const& rOp : WavefunctionDM)
   {
      for (auto const& cOp : rOp)
      {
         int tp;
         LinearAlgebra::Range rtp;
         std::tie(tp, rtp) = B.Lookup(rOp.row()());
         int t;
         LinearAlgebra::Range rt;
         std::tie(t, rt) = B.Lookup(cOp.col()());
         CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal")(B[tp])(B[t])(WavefunctionDM)(B);
         PsiDMList[tp](rtp, rt) = *J;
      }
   }

   // diagonalize the density matrix
   this->DiagonalizeDMHelper(false);

   // Get the eigenvalues with respect to WavefunctionDM
   std::size_t i = 0;
   for (std::size_t j = 0; j < PsiDMList.size(); ++j)
   {
      RawDMType M = RawDMList[j] * PsiDMList[j] * herm(RawDMList[j]);
      for (std::size_t k = 0; k < M.size1(); ++k)
      {
         EigenInfoList[i].Eigenvalue = LinearAlgebra::real(M(k,k));
         //         TRACE(EigenInfoList[i].Eigenvalue);
         ++i;
      }
   }
   CHECK_EQUAL(i, EigenInfoList.size());
   std::sort(EigenInfoList.begin(), EigenInfoList.end(), EigenCompare);
}
#endif

//
// DensityMatrix<SimpleOperator>
//

DensityMatrix<SimpleOperator>::DensityMatrix(SimpleOperator const& Op)
   : B(Op.Basis1())
{
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), Op.Basis2());
   // initialize the RawDMList.  At the same time, get the maximum linear size of the subspaces
   RawDMList.clear();
   RawDMList.reserve(B.size());
   MaxLinearDimension = 0;
   for (std::size_t q = 0; q < B.size(); ++q)
   {
      int Dim = B.dim(q);
      //      std::cout << "dimension of " << q << " is " << Dim << std::endl;
      RawDMList.emplace_back(Dim, Dim, 0);
      MaxLinearDimension = std::max<int>(MaxLinearDimension, Dim);
   }
   // Fill the raw density matrices
   for (auto const& rOp : Op)
   {
      for (auto const& cOp : rOp)
      {
         int tp, rtp;
         std::tie(tp, rtp) = B.Lookup(rOp.row());
         int t, rt;
         std::tie(t, rt) = B.Lookup(cOp.col());
         CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal");
         RawDMList[tp](rtp, rt).set_wait(cOp.value);
      }
   }
   // diagonalize them
   this->DiagonalizeDMHelper();
}

void UpdateKeepList(KeepListType& KeepList,
                    std::set<QuantumNumbers::QuantumNumber> const& SiteQN,
                    VectorBasis const& FullBasis,
                    std::list<EigenInfo>& KeepStates,
                    std::list<EigenInfo>& DiscardStates,
                    TruncationInfo& Info)
{
   typedef std::set<QuantumNumbers::QuantumNumber> SiteQNType;
   typedef std::list<EigenInfo> StatesListType;

   // Make a set of all quantum numbers that are already kept
   typedef std::set<QuantumNumbers::QuantumNumber> qnType;
   qnType KeptQN;
   for (StatesListType::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
      KeptQN.insert(FullBasis[I->Subspace]);

   // Iterate through the KeepList and get the transformed states, make sure they
   // exist in KeepStates
   for (KeepListType::const_iterator I = KeepList.begin(); I != KeepList.end(); ++I)
   {
      qnType qn;  // the set of updated quantum numbers corresponding to keep state I
      for (SiteQNType::const_iterator Q = SiteQN.begin(); Q != SiteQN.end(); ++Q)
         transform_targets(*I, *Q, std::inserter(qn, qn.end()));
      // make sure that at least one of qn exists in KeptQN
      SiteQNType Intersect;
      std::set_intersection(KeptQN.begin(), KeptQN.end(), qn.begin(), qn.end(),
                            std::inserter(Intersect, Intersect.end()));
      if (Intersect.empty())
      {
         // No kept state, go through the discard list and resurrect a state
         StatesListType::iterator Piv = DiscardStates.begin();
         while (Piv != DiscardStates.end() && qn.count(FullBasis[Piv->Subspace]) == 0)
            ++Piv;
         if (Piv == DiscardStates.end())
         {
            WARNING("Unexpected: cannot force quantum number into basis")(*I);
         }
         else
         {
            DEBUG_WARNING("Forcing keep state")(*I)(FullBasis[Piv->Subspace])(Piv->Weight());
            ++Info.KeptStates_;
            Info.KeptWeight_ += Piv->Weight();
            Info.KeptEntropy_ += Piv->Entropy(Info.TotalWeight_);
            KeptQN.insert(FullBasis[Piv->Subspace]);
            KeepStates.push_back(*Piv);
            DiscardStates.erase(Piv);
         }
      }
   }

   // Now construct the new KeepList.
   KeepList.clear();
   for (StatesListType::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
   {
      KeepList.insert(FullBasis[I->Subspace]);
   }
   //std::copy(KeepList.begin(),KeepList.end(),std::ostream_iterator<QuantumNumber>(std::cout, " "));
   //std::cout << '\n';
}

//
// Singular value decomposition
//

SingularDecompositionBase::SingularDecompositionBase()
{
}

SingularDecompositionBase::~SingularDecompositionBase()
{
}

void SingularDecompositionBase:: Diagonalize(std::vector<RawDMType>&& M)
{
   ESum = 0;  // Running sum of squares of the singular values
   for (std::size_t i = 0; i < M.size(); ++i)
   {
      int nv = std::min(M[i].rows(), M[i].cols());
      Matrix U(M[i].rows(),nv);
      Matrix Vh(nv, M[i].cols());
      RealVector D_device(nv);
      SVD(std::move(M[i]), U, D_device, Vh);
      cpu::RealVector D = get_wait(D_device);

      LeftVectors.push_back(std::move(U));
      RightVectors.push_back(std::move(Vh));
      SingularValues.push_back(std::move(D_device));

      int CurrentDegree = degree(this->Lookup(i));
      for (unsigned j = 0; j < D.size(); ++j)
      {
         double Weight = D[j]*D[j];
         EigenInfoList.push_back(EigenInfo(Weight, CurrentDegree, i, j));
         ESum += Weight * CurrentDegree;
      }
   }

   std::sort(EigenInfoList.begin(), EigenInfoList.end(), EigenCompare);
}


SingularDecomposition<MatrixOperator, MatrixOperator>::SingularDecomposition(MatrixOperator const& M)
   : B1(M.Basis1()), B2(M.Basis2()), IndexOfi(B1.size(), -1)
{
   // Assemble the raw matrices

   // Iterate through Basis1 and Basis2 and assemble the list of used quantum numbers, in order,
   // and how they map onto the components in Basis1 and Basis2.
   // We can make use of the fact that each quantum number only appears at most once in
   // the LinearBasis.
   for (unsigned i = 0; i < B1.size(); ++i)
   {
      QuantumNumber q = B1[i];
      int j = B2.find_first(q);
      if (j == -1)
         continue;

      IndexOfi[i] = UsedQuantumNumbers.size();
      q_iLinear.push_back(i);
      q_jLinear.push_back(j);
      UsedQuantumNumbers.push_back(q);
   }

   // Now assemble the list of raw matrices, initialized to zero
   std::vector<RawDMType> Matrices;
   Matrices.reserve(q_iLinear.size());
   for (unsigned n = 0; n < q_iLinear.size(); ++n)
   {
      Matrices.emplace_back(B1.dim(q_iLinear[n]), B2.dim(q_jLinear[n]), 0.0);
   }

   // fill the raw matrices with the components in M
   for (auto const& Mr : M)
   {
      std::pair<int, blas::Range> iIndex = B1.Lookup(Mr.row());
      for (auto const& Mc : Mr)
      {
                  // determine where this (i,j) component fits within the matrices
         std::pair<int, blas::Range> jIndex = B2.Lookup(Mc.col());
         // map the index into the used subspaces
         int Subspace = IndexOfi[iIndex.first];
         CHECK(Subspace != -1);
         CHECK_EQUAL(q_jLinear[Subspace], jIndex.first);
         Matrices[Subspace](iIndex.second, jIndex.second) = Mc.value;
      }
   }

   // do the SVD
   this->Diagonalize(std::move(Matrices));
}

QuantumNumber
SingularDecomposition<MatrixOperator, MatrixOperator>::Lookup(int Subspace) const
{
   return UsedQuantumNumbers[Subspace];
}

std::tuple<MatrixOperator, RealDiagonalOperator, MatrixOperator>
SingularDecomposition<MatrixOperator, MatrixOperator>::
ConstructOrthoMatrices(std::vector<std::vector<int> > const& LinearMapping)
{
   // Construct the truncated basis
   int NumQ = UsedQuantumNumbers.size();
   VectorBasis NewBasis(B1.GetSymmetryList());
   std::vector<int> NewSubspace(NumQ, -1); // maps the q to the new label
   for (int q = 0; q < NumQ; ++q)
   {
      if (LinearMapping[q].size() > 0)
      {
         NewSubspace[q] = NewBasis.size();
         NewBasis.push_back(UsedQuantumNumbers[q], LinearMapping[q].size());
      }
   }

   // Now we construct the actual matrices
   MatrixOperator A = MatrixOperator(B1.MappedBasis(), NewBasis);
   MatrixOperator B = MatrixOperator(NewBasis, B2.MappedBasis());
   for (int q = 0; q < NumQ; ++q)
   {
      int ss = NewSubspace[q];
      if (ss == -1) // is this subspace used?
         continue;

      QuantumNumber Q = UsedQuantumNumbers[q];

      // The set of indices of states we need to keep
      std::vector<int> lm(LinearMapping[q].begin(), LinearMapping[q].end());

      // Now assemble this quantum number component for A
      int i = B1.MappedBasis().find_first(Q);
      while (i != -1)
      {
         Matrix Temp(B1.Lookup(i).second.size(), lm.size());
         assign_slice(Temp, LeftVectors[q], B1.Lookup(i).second, lm);
         A.insert(i,ss, std::move(Temp));
         i = B1.MappedBasis().find_next(Q, i);
      }

      // and B
      int j = B2.MappedBasis().find_first(Q);
      while (j != -1)
      {
         Matrix Temp(lm.size(), B2.Lookup(j).second.size());
         assign_slice(Temp, RightVectors[q], lm, B2.Lookup(j).second);
         B.insert(ss,j, std::move(Temp));
         j = B2.MappedBasis().find_next(Q, j);
      }
   }

   // Finally the center matrix
   RealDiagonalOperator C = RealDiagonalOperator(NewBasis, NewBasis);
   for (int i = 0; i < NumQ; ++i)
   {
      if (!LinearMapping[i].empty())
      {
         int b = NewSubspace[i];
         RealDiagonalMatrix Temp(int(LinearMapping[i].size()));
         assign_permutation(Temp.diagonal(), SingularValues[i], &LinearMapping[i][0]);
         C.insert(b,b, std::move(Temp));
      }
   }

   A.debug_check_structure();
   B.debug_check_structure();
}

QuantumNumber
SingularDecomposition<StateComponent, StateComponent>::Lookup(int Subspace) const
{
   return UsedQuantumNumbers[Subspace];
}

SingularDecomposition<StateComponent, StateComponent>::
SingularDecomposition(StateComponent const& A, ProductBasis<BasisList, BasisList> Factors_)
   : SingularDecompositionBase(), B1(A.Basis1()), B2(A.Basis2()), Factors(std::move(Factors_))
{
   using QuantumNumbers::QuantumNumberList;

   CHECK_EQUAL(A.LocalBasis(), Factors.Basis());

   // denote the operator by <j' || A[k] || j>
   // where A[k] decomposes as S[k_1] T[k_2]
   // The matrix is block diagonal with respect to (j' - k_1) and (j + k_2).
   // So the quantum number for the singular values runs over these numbers.

   // map each distinct quantum number onto successive integers
   int NumQuantum = 0;
   std::map<QuantumNumber, int> QuantumNumberLinearMapping;// maps to the index in UsedQuantumNumbers

   // linear dimensions of the subspaces for quantum number pair (j' - k_1) = (j + k_2) = r
   std::vector<std::pair<int, int> > LinearDimensions;

   // for Basis1, construct the combined basis
   for (unsigned k1 = 0; k1 < Factors.Left().size(); ++k1)
   {
      QuantumNumber k1q_bar = adjoint(Factors.Left()[k1]);

      for (unsigned i = 0; i < A.Basis1().size(); ++i)
      {
         QuantumNumberList ql = transform_targets(A.Basis1()[i], k1q_bar);
         for (QuantumNumberList::const_iterator qi = ql.begin(); qi != ql.end(); ++qi)
         {
            int LinearIndex;
            std::map<QuantumNumber, int>::iterator I = QuantumNumberLinearMapping.find(*qi);
            if (I == QuantumNumberLinearMapping.end())
            {
               // not found, add it
               LinearIndex = UsedQuantumNumbers.size();
               QuantumNumberLinearMapping[*qi] = LinearIndex;
               UsedQuantumNumbers.push_back(*qi);
               LinearDimensions.push_back(std::pair<int, int>(0,0));
            }
            else
               LinearIndex = I->second;

            int dim = A.Basis1().dim(i);
            int current = LinearDimensions[LinearIndex].first;
            LeftSubspaceInfo.push_back(ProductSubspaceInfo(k1, i, LinearIndex,
                                                           blas::Range(current, current+dim)));
            LinearDimensions[LinearIndex].first += dim;
         }
      }
   }

   // same for Basis2

   for (unsigned k2 = 0; k2 < Factors.Right().size(); ++k2)
   {
      QuantumNumber k2q = Factors.Right()[k2];

      for (unsigned i = 0; i < A.Basis2().size(); ++i)
      {
         QuantumNumberList ql = transform_targets(A.Basis2()[i], k2q);

         for (QuantumNumberList::const_iterator qi = ql.begin(); qi != ql.end(); ++qi)
         {
            int LinearIndex;
            std::map<QuantumNumber, int>::iterator I = QuantumNumberLinearMapping.find(*qi);
            if (I == QuantumNumberLinearMapping.end())
            {
               // not found, add it.  This is a degenerate case in this context,
               // since the quantum number doesn't exist in the left basis so
               // this component is a zero-dimensional matrix
               LinearIndex = UsedQuantumNumbers.size();
               QuantumNumberLinearMapping[*qi] = LinearIndex;
               UsedQuantumNumbers.push_back(*qi);
               LinearDimensions.push_back(std::pair<int, int>(0,0));
            }
            else
               LinearIndex = I->second;

            int dim = A.Basis2().dim(i);
            int current = LinearDimensions[LinearIndex].second;
            RightSubspaceInfo.push_back(ProductSubspaceInfo(k2, i, LinearIndex,
                                                            blas::Range(current, current+dim)));
            LinearDimensions[LinearIndex].second += dim;
         }
      }
   }

   NumQuantum = UsedQuantumNumbers.size();

   // Now we can construct the actual matrices, initialized to zero
   std::vector<RawDMType> Matrices;
   Matrices.reserve(NumQuantum);
   for (int i = 0; i < NumQuantum; ++i)
      Matrices.emplace_back(LinearDimensions[i].first, LinearDimensions[i].second, 0.0);

   // and fill them
   for (unsigned lin_1 = 0; lin_1 < LeftSubspaceInfo.size(); ++lin_1)
   {
      for (unsigned lin_2 = 0; lin_2 < RightSubspaceInfo.size(); ++lin_2)
      {
         // we only care about blocks that are on the diagonal
         if (LeftSubspaceInfo[lin_1].LinearIndex != RightSubspaceInfo[lin_2].LinearIndex)
            continue;

         // and make sure that block also has some non-zero element
         if (LeftSubspaceInfo[lin_1].LinearRange.size() == 0 || RightSubspaceInfo[lin_2].LinearRange.size() == 0)
            continue;

         // sum over k
         ProductBasis<BasisList, BasisList>::const_iterator klEnd
            = Factors.end(LeftSubspaceInfo[lin_1].k,
                          RightSubspaceInfo[lin_2].k);
         ProductBasis<BasisList, BasisList>::const_iterator klIter
            = Factors.begin(LeftSubspaceInfo[lin_1].k,
                            RightSubspaceInfo[lin_2].k);
         for ( ; klIter != klEnd; ++klIter)
         {
            auto J = A[*klIter].data().row(LeftSubspaceInfo[lin_1].s).find(RightSubspaceInfo[lin_2].s);
            if (J != A[*klIter].data().row(LeftSubspaceInfo[lin_1].s).end())
            {
               double Coeff
                  = inverse_product_coefficient(Factors.Left()[LeftSubspaceInfo[lin_1].k],
                                                Factors.Right()[RightSubspaceInfo[lin_2].k],
                                                Factors.Basis()[*klIter],
                                                A.Basis1()[LeftSubspaceInfo[lin_1].s],
                                                A.Basis2()[RightSubspaceInfo[lin_2].s],
                                                UsedQuantumNumbers[LeftSubspaceInfo[lin_1].LinearIndex]);

               // The normalization factor is to satisfy the left-normalization constraint
               // with the quantum-number dependent prefactor.  See section 7.1 of the nonabelianmp,
               // equation eq:LeftAMatrixNorm.  Given an ordinary unitary matrix, if we insert it
               // into an A-matrix, the normalization changes due to the conventions of the reduced
               // matrix elements.  To counteract this, we need to scale the unitary matrix by
               // a factor as it is inserted into the A-matrix, and apply the inverse of the
               // scale factor at this point here.
               double NormFactor = std::sqrt(double(degree(A.Basis1()[LeftSubspaceInfo[lin_1].s]))
                                       / degree(UsedQuantumNumbers[LeftSubspaceInfo[lin_1].LinearIndex]));

               Matrices[LeftSubspaceInfo[lin_1].LinearIndex](LeftSubspaceInfo[lin_1].LinearRange,
                                                             RightSubspaceInfo[lin_2].LinearRange)
                  += NormFactor * Coeff * (*J).value;
            }
         }
      }
   }

   // do the SVD
   this->Diagonalize(std::move(Matrices));
}

std::tuple<StateComponent, RealDiagonalOperator, StateComponent>
SingularDecomposition<StateComponent, StateComponent>::
ConstructOrthoMatrices(std::vector<std::vector<int> > const& LinearMapping)
{
   int NumQ = UsedQuantumNumbers.size();
   VectorBasis NewBasis(B1.GetSymmetryList());

   // Now we construct the truncated basis
   std::vector<int> NewSubspace(NumQ, -1); // maps the q to the new label
   for (int q = 0; q < NumQ; ++q)
   {
      if (LinearMapping[q].size() > 0)
      {
         NewSubspace[q] = NewBasis.size();
         NewBasis.push_back(this->Lookup(q), LinearMapping[q].size());
      }
   }
   // Now we construct the actual matrices, firstly for the left
   StateComponent A(Factors.Left(), B1, NewBasis);
   for (unsigned i = 0; i < LeftSubspaceInfo.size(); ++i)
   {
      int ss = LeftSubspaceInfo[i].LinearIndex;
      if (!LinearMapping[ss].empty())
      {
         std::vector<int> lm(LinearMapping[ss].begin(), LinearMapping[ss].end());

         // Apply the normalization factor to ensure A satisfies the left-orthonormalization constraint.
         // This is the reciprocal of the factor we applied previously.
         double NormFactor = std::sqrt(double(degree(A.Basis2()[NewSubspace[ss]]))
                                       / degree(A.Basis1()[LeftSubspaceInfo[i].s]));

         Matrix Temp(LeftSubspaceInfo[i].LinearRange.size(), lm.size());
         assign_slice(Temp,  LeftVectors[ss], LeftSubspaceInfo[i].LinearRange, lm);
         scale(Temp, complex(NormFactor));
         A[LeftSubspaceInfo[i].k].insert(LeftSubspaceInfo[i].s, NewSubspace[ss], std::move(Temp));
      }
   }

   // now the right
   StateComponent B(Factors.Right(), NewBasis, B2);
   for (unsigned i = 0; i < RightSubspaceInfo.size(); ++i)
   {
      int ss = RightSubspaceInfo[i].LinearIndex;
      std::vector<int> lm(LinearMapping[ss].begin(), LinearMapping[ss].end());
      if (!lm.empty())
      {
         Matrix Temp(lm.size(), RightSubspaceInfo[i].LinearRange.size());
         assign_slice(Temp, RightVectors[ss], lm, RightSubspaceInfo[i].LinearRange);
         B[RightSubspaceInfo[i].k].insert(NewSubspace[ss], RightSubspaceInfo[i].s, Temp);
      }
   }

   // Finally the center matrix
   RealDiagonalOperator C(NewBasis, NewBasis);
   for (int i = 0; i < NumQ; ++i)
   {
      if (!LinearMapping[i].empty())
      {
         int b = NewSubspace[i];
         RealDiagonalMatrix Temp(int(LinearMapping[i].size()));
         assign_permutation(Temp.diagonal(), SingularValues[i], &LinearMapping[i][0]);
         C.insert(b,b, std::move(Temp));
      }
   }
}
