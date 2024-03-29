// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mps/density.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(DENSITY_H_FDSHYFUIH38348UER8J8U3)
#define DENSITY_H_FDSHYFUIH38348UER8J8U3

#include "mps/state_component.h"
#include "truncation.h"
#include "linearalgebra/diagonalmatrix.h"
#include <tuple>
#include <set>
#include <list>

// sometimes we need a number of states that stands for 'infinite'
int const DefaultMaxStates = 50000;

typedef std::set<QuantumNumbers::QuantumNumber> KeepListType;

// For the quantum numbers in KeepList, make sure that at least one state
// in KeepList \otimes SiteQN is represented in the KeepStates.
// States are shifted from DiscardedStates to KeepStates in order to satisfy
// this condition.
// FullBasis is the actual basis (kept + discarded states), which is used to
// get the quantum number from the subspace index of the kept and discarded states.
// The updated KeepList is then set to the quantum numbers in the final basis.
//
// The idea behind the KeepList is that quantum number sectors are reachable.
// For example, suppose we are doing a left-to-right sweep and updating
// some site A^s.  The KeepList will be the set of quantum number sectors
// in A.Basis1() that we want to keep (i.e. they are part of the kept part
// of the wavefunction.  This excludes states that might be in the basis
// for other reasons, eg subspace expansion).  Now when we truncate Basis2()
// of A^s, we want to ensure that Basis2() * s covers every quantum number
// sector in KeepList.  This guarantees that on the subsequent right-to-left
// sweep, all of the kept sectors have a possibility to be included in the
// wavefunction.
void UpdateKeepList(KeepListType& KeepList,
                    std::set<QuantumNumbers::QuantumNumber> const& SiteQN,
                    VectorBasis const& FullBasis,
                    std::list<EigenInfo>& KeepStates,
                    std::list<EigenInfo>& DiscardStates,
                    TruncationInfo& Info);

// Function to get the dimension of each quantum number subspace in the
// eigenvalue list

typedef std::map<QuantumNumbers::QuantumNumber, int> EigenDimensionsType;

template <typename FwdIter>
EigenDimensionsType
EigenDimensions(VectorBasis const& B, FwdIter first, FwdIter last)
{
   std::map<QuantumNumbers::QuantumNumber, int> Result;
   while (first != last)
   {
      ++Result[B[first->Subspace]];
      ++first;
   }
   return Result;
}

// Adjusts the dimension of each quantum number to be 1+floor((Original-1)*Fraction+0.5)
inline
void UpdateEigenDimensions(EigenDimensionsType& EDim,
                           double Fraction)
{
   for (EigenDimensionsType::iterator I = EDim.begin(); I != EDim.end(); ++I)
      I->second = 1 + std::max(0, int(floor((I->second-1)*Fraction + 0.5)));
}

// Determine the states to keep from the list of dimensions
template <typename FwdIter>
TruncationInfo
TruncateFixEigenDimensions(EigenDimensionsType const& EDim, FwdIter first, FwdIter last);

template <typename FwdIter>
std::list<EigenInfo>
TruncateFixEigenDimensions(EigenDimensionsType const& EDim_, VectorBasis const& B,
                           int MinStates,
                           FwdIter first, FwdIter last,
                           TruncationInfo& Info)
{
   std::list<EigenInfo> Result;
   EigenDimensionsType EDim = EDim_;  // make a copy, so we can modify it locally
   Info.TotalWeight_ = DensitySigma(first, last);
   Info.TotalStates_ = std::distance(first, last);
   Info.TotalEntropy_ = DensityEntropy(first, last, Info.TotalWeight_);
   Info.KeptStates_ = 0;
   Info.KeptWeight_ = 0;
   Info.KeptEntropy_ = 0;
   while (first != last)
   {
      if (--EDim[B[first->Subspace]] >= 0 || MinStates > 0)
      {
         Result.push_back(*first);
         ++Info.KeptStates_;
         Info.KeptWeight_ += first->Weight();
         Info.KeptEntropy_ += first->Entropy(Info.TotalWeight_);
         --MinStates;
      }
      ++first;
   }
   return Result;
}

//
// LinearBasis
//
// The LinearBasis class is used to map a discontiguous basis
// into a contiguous one, where each quantum number occurs only once.
// This is similar in principle to the regularize() operation on
// irred tensors - indeed, regularize() could (and possibly should)
// be implemented via a LinearBasis.
// For each index i into the original basis, Lookup(i) returns a pair,
// being the new subspace number (that uniquely identifies a symmetry sector),
// and the location within the subspace where state i belongs.  For a VectorBasis,
// this is a LinearAlgebra::Range.  For a BasisList, this is a simple integer.
//

template <typename BasisT>
class LinearBasis;

template <>
class LinearBasis<VectorBasis> : public VectorBasis
{
   public:
      typedef VectorBasis MappedBasisType;
      typedef std::pair<int, LinearAlgebra::Range> MapType;

      explicit LinearBasis(MappedBasisType const& B);

      MappedBasisType const& MappedBasis() const { return Orig; }

      MapType Lookup(int Original) const { return Mapping[Original]; }

   private:
      MappedBasisType Orig;
      LinearAlgebra::Vector<std::pair<int, LinearAlgebra::Range> > Mapping;
};

template <>
class LinearBasis<BasisList> : public VectorBasis
{
   public:
      typedef BasisList MappedBasisType;
      typedef std::pair<int, int> MapType;

      explicit LinearBasis(MappedBasisType const& B);

      MappedBasisType const& MappedBasis() const { return Orig; }

      MapType Lookup(int Original) const { return Mapping[Original]; }

      int ReverseLookup(int s, int index) const;

   private:
      MappedBasisType Orig;
      LinearAlgebra::Vector<std::pair<int, int> > Mapping;
};

class DensityMatrixBase
{
   public:
      typedef std::vector<EigenInfo> EigenInfoListType;
      typedef EigenInfoListType::const_iterator const_iterator;

      // iterators for the beginning and end of the sorted list of eigenvalues.
      const_iterator begin() const { return EigenInfoList.begin(); }
      const_iterator end() const { return EigenInfoList.end(); }

      int size() const { return EigenInfoList.size(); }

      // shows a report of the density matrix eigenvalues, showing at most MaxEigenvalues of
      // the eigenvalue, cumulative truncation error,
      // TODO: fix the code rot
      // Base2 means show entropy as base 2 rather than natural log
      // ShowDegen shows multiplets as repeated eigenvalues
      std::ostream& DensityMatrixReport(std::ostream& out, int MaxEigenvalues = -1, bool Base2 = false,
					bool ShowDegen = false, bool Quiet = false);

      // returns the sum of the eigenvalues
      double EigenSum() const { return ESum; }

      double Entropy(bool Base2 = false) const;

      double EvaluateCasimir(int n) const;

      double EvaluateCasimirMoment(int n) const;

   protected:
      DensityMatrixBase() {}

      // this function diagonalizes the RawDMList, and initializes the EigenInfoList.
      // It is assumed
      // on entry that Basis is valid, and RawDMList is the undiagonalized density matrix.
      void DiagonalizeDMHelper(bool Sort = true);

      virtual ~DensityMatrixBase() {}

   public:
      virtual QuantumNumber Lookup(int Subspace) const = 0;

   protected:
      // RawDM is (new, old)
      typedef LinearAlgebra::Matrix<std::complex<double> > RawDMType;

      std::vector<RawDMType> RawDMList;
      std::vector<EigenInfo> EigenInfoList;
      int MaxLinearDimension;
      double ESum;
};

template <typename OperatorT>
class DensityMatrix;

template <>
class DensityMatrix<MatrixOperator> : public DensityMatrixBase
{
   public:
      typedef MatrixOperator OperatorType;
      typedef OperatorType::basis1_type BasisType;

      // constructs the density matrix eigenstates from Op.  Op must be symmetric (Hermitian)
      DensityMatrix(OperatorType const& Op);

      // constructs the density matrix eigenstates from Op, but uses WavefunctionDM
      // to get the actual eigenvalues.  This is for cases where the density matrix
      // inludes a mixing factor.
      DensityMatrix(OperatorType const& Op, OperatorType const& WavefunctionDM);

      LinearBasis<BasisType> const& Basis() const { return B; }

      // constructs a truncation operator that projects onto the given eigenstates.
      // The resulting operator has Basis1()' = truncated basis, Basis2()' = original basis
      template <typename FwdIter>
      OperatorType ConstructTruncator(FwdIter First, FwdIter Last) const;

   private:

   public:
      QuantumNumber Lookup(int Subspace) const { return B[Subspace]; }

   private:
      LinearBasis<BasisType> B;
};

template <>
class DensityMatrix<SimpleOperator> : public DensityMatrixBase
{
   public:
      typedef SimpleOperator OperatorType;
      typedef OperatorType::basis1_type BasisType;

      // constructs the density matrix eigenstates from Op.  Op must be symmetric (Hermitian)
      DensityMatrix(OperatorType const& Op);

      LinearBasis<BasisType> const& Basis() const { return B; }

      // constructs a truncation operator that projects onto the given eigenstates
      template <typename FwdIter>
      OperatorType ConstructTruncator(FwdIter First, FwdIter Last) const;

      template <typename FwdIter>
      OperatorType ConstructUnnormalizedTruncator(FwdIter First, FwdIter Last) const;

   private:
      QuantumNumber Lookup(int Subspace) const { return B[Subspace]; }

      LinearBasis<BasisType> B;
};

//
// Singular value decomposition
//
// These templates have two parameters, being the types of the U and V^\dagger matrices
// that we want to construct.  The eigenvalues in the EigenInfo list are squared so that
// they properly represent weights.  This means we can do a drop-in replacement of
// the density matrix, where appropriate.
//

template <typename Mat1, typename Mat2>
class SingularDecomposition;

class SingularDecompositionBase
{
   public:
      typedef std::vector<EigenInfo> EigenInfoListType;
      typedef EigenInfoListType::const_iterator const_iterator;

      // iterators for the beginning and end of the sorted list of eigenvalues.
      const_iterator begin() const { return EigenInfoList.begin(); }
      const_iterator end() const { return EigenInfoList.end(); }

      double EigenSum() const { return ESum; }

      enum WhichVectors { Left, Right, Both };

      std::ostream& DensityMatrixReport(std::ostream& out);

   protected:
      typedef LinearAlgebra::Matrix<std::complex<double> > RawDMType;

      SingularDecompositionBase();

      void Diagonalize(std::vector<RawDMType> const& M, WhichVectors Which = Both, WhichVectors WhichCalculate = Both);


   public:
      virtual QuantumNumber Lookup(int Subspace) const = 0;

   protected:
      virtual ~SingularDecompositionBase();

      std::vector<RawDMType> LeftVectors, RightVectors;
      std::vector<EigenInfo> EigenInfoList;
      std::vector<LinearAlgebra::Vector<double>> SingularValues;  // to avoid taking sqrt of density eigenvalues
      int MaxLinearDimension;
      double ESum;
};

struct ProductSubspaceInfo
{
   int k;  // the local basis index
   int s;  // index into the matrix basis

   int LinearIndex;                  // index into the combined linear basis
   LinearAlgebra::Range LinearRange; // corresponding range of the linear basis

   ProductSubspaceInfo() {}
   ProductSubspaceInfo(int k_, int s_, int LinearIndex_, LinearAlgebra::Range LinearRange_)
      : k(k_), s(s_), LinearIndex(LinearIndex_), LinearRange(LinearRange_) {}
};

template <>
class SingularDecomposition<MatrixOperator, MatrixOperator> : public SingularDecompositionBase
{
   public:
      typedef MatrixOperator left_type;
      typedef MatrixOperator right_type;
      typedef RealDiagonalOperator diagonal_type;

      explicit SingularDecomposition(MatrixOperator const& M);

      // Optionally, if we want to construct the full basis for the left or right, we can do that
      SingularDecomposition(MatrixOperator const& M, WhichVectors Which);

      // Optionally, also choose which singular vectors to calculate
      SingularDecomposition(MatrixOperator const& M, WhichVectors Which, WhichVectors WhichCalculate);

      template <typename FwdIter>
      void ConstructMatrices(FwdIter first, FwdIter last, MatrixOperator& A, RealDiagonalOperator& C, MatrixOperator& B);

      // Construct only the left singular vectors
      template <typename FwdIter>
      MatrixOperator
      ConstructLeftVectors(FwdIter first, FwdIter last);

      // Construct only the right singular vectors
      template <typename FwdIter>
      MatrixOperator
      ConstructRightVectors(FwdIter first, FwdIter last);

      // Construct the diagonal matrix of singular values
      template <typename FwdIter>
      RealDiagonalOperator
      ConstructSingularValues(FwdIter first, FwdIter last);


      LinearBasis<VectorBasis> Basis1() const { return B1; }
      LinearBasis<VectorBasis> Basis2() const { return B2; }

   private:
      template <typename FwdIter>
      std::tuple<VectorBasis, std::vector<std::set<int>>, std::vector<int>>
      ConstructMapping(FwdIter first, FwdIter last);

      MatrixOperator
      DoConstructLeftVectors(std::tuple<VectorBasis, std::vector<std::set<int>>, std::vector<int>> const Mapping);

      MatrixOperator
      DoConstructRightVectors(std::tuple<VectorBasis, std::vector<std::set<int>>, std::vector<int>> const Mapping);

      RealDiagonalOperator
      DoConstructSingularValues(std::tuple<VectorBasis, std::vector<std::set<int>>, std::vector<int>> const Mapping);

      QuantumNumber Lookup(int Subspace) const;

      LinearBasis<VectorBasis> B1, B2;
      std::vector<int> q_iLinear, q_jLinear;
      std::vector<int> IndexOfi;
      std::vector<QuantumNumber> UsedQuantumNumbers;  // in order

      void ConstructOrthoMatrices(std::vector<std::set<int> > const& LinearMapping,
                                  MatrixOperator& A, RealDiagonalOperator& C, MatrixOperator& B);
};

typedef SingularDecomposition<MatrixOperator, MatrixOperator> CMatSVD;  // avoid typing...

template <>
class SingularDecomposition<StateComponent, StateComponent> : public SingularDecompositionBase
{
   public:
      typedef StateComponent left_type;
      typedef StateComponent right_type;
      typedef RealDiagonalOperator diagonal_type;

      SingularDecomposition(StateComponent const& A, ProductBasis<BasisList, BasisList> const& Factors);

      template <typename FwdIter>
      void ConstructMatrices(FwdIter first, FwdIter last,
                             StateComponent& A,
                             RealDiagonalOperator& C,
                             StateComponent& B);

   private:
      QuantumNumber Lookup(int Subspace) const;

      void ConstructOrthoMatrices(std::vector<std::set<int> > const& LinearMapping,
                                  StateComponent& A, RealDiagonalOperator& C, StateComponent& B);

      VectorBasis B1, B2;
      ProductBasis<BasisList, BasisList> Factors;
      std::vector<ProductSubspaceInfo> LeftSubspaceInfo, RightSubspaceInfo;
      std::vector<QuantumNumber> UsedQuantumNumbers;  // in order

   //      LinearBasis<VectorBasis> Basis1, Basis2;
};

typedef SingularDecomposition<StateComponent, StateComponent> AMatSVD;  // avoid typing...

#include "density.cc"

#endif
