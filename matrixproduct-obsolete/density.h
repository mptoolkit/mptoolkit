// -*- C++ -*- $Id$

#if !defined(DENSITY_H_FDSHYFUIH38348UER8J8U3)
#define DENSITY_H_FDSHYFUIH38348UER8J8U3

#include "mps/state_component.h"
#include "linearalgebra/diagonalmatrix.h"
#include <boost/tuple/tuple.hpp>
#include <set>
#include <list>

// sometimes we need a number of states that stands for 'infinite'
int const DefaultMaxStates = 50000;  

typedef IrredTensor
        <
           LinearAlgebra::DiagonalMatrix<double>
        ,  VectorBasis
        ,  VectorBasis
        >
        RealDiagonalOperator;

typedef std::set<QuantumNumbers::QuantumNumber> KeepListType;

// For the quantum numbers in KeepList, make sure that at least one state
// in KeepList \otimes SiteQN is represented in the KeepStates.
// The updated KeepList is then set to the quantum numbers in the final basis.
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

      // shows a report of the density matrix eigenvalues, showing at most MaxEigenvalues of
      // the eigenvalue, cumulative truncation error, 
      // TODO: fix the code rot
      std::ostream& DensityMatrixReport(std::ostream& out, int MaxEigenvalues = -1, bool Base2 = false);

      // returns the sum of the eigenvalues
      double EigenSum() const { return ESum; }

      double EvaluateCasimir(int n) const;

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

      // constructs a truncation operator that projects onto the given eigenstates
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

   protected:
      typedef LinearAlgebra::Matrix<std::complex<double> > RawDMType;

      SingularDecompositionBase();

      void Diagonalize(std::vector<RawDMType> const& M);

   public:
      virtual QuantumNumber Lookup(int Subspace) const = 0;

   protected:
      virtual ~SingularDecompositionBase();

      std::vector<RawDMType> LeftVectors, RightVectors;
      std::vector<EigenInfo> EigenInfoList;
      std::vector<LinearAlgebra::Vector<double> > SingularValues;  // to avoid taking sqrt of density eigenvalues
      int MaxLinearDimension;
      double ESum;
};

template <>
class SingularDecomposition<MPStateComponent, MatrixOperator> : public SingularDecompositionBase
{
   public:
      typedef MPStateComponent left_type;
      typedef MatrixOperator right_type;
      typedef DiagonalOperator diagonal_type;

      SingularDecomposition(MPStateComponent const& A);

#if 0
      // constructs the U, V^\dagger matrices
      template <typename FwdIter>
      boost::tuple<left_type, diagonal_type, right_type> ConstructTruncator(FwdIter First, FwdIter Last) const;
#endif

   private:
      LinearBasis<VectorBasis> Basis1, Basis2;
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
class SingularDecomposition<MPStateComponent, MPStateComponent> : public SingularDecompositionBase
{
   public:
      typedef MPStateComponent left_type;
      typedef MPStateComponent right_type;
      typedef RealDiagonalOperator diagonal_type;

      SingularDecomposition(MPStateComponent const& A, ProductBasis<BasisList, BasisList> const& Factors);

      template <typename FwdIter>
      void ConstructMatrices(FwdIter first, FwdIter last, 
			     MPStateComponent& A, 
			     RealDiagonalOperator& C, 
			     MPStateComponent& B);

   private:
      QuantumNumber Lookup(int Subspace) const;
   
      void ConstructOrthoMatrices(std::vector<std::set<int> > const& LinearMapping,
				  MPStateComponent& A, MatrixOperator& C, MPStateComponent& B);

      VectorBasis B1, B2;
      ProductBasis<BasisList, BasisList> Factors;
      std::vector<ProductSubspaceInfo> LeftSubspaceInfo, RightSubspaceInfo;
      std::vector<QuantumNumber> UsedQuantumNumbers;  // in order

   //      LinearBasis<VectorBasis> Basis1, Basis2;
};

typedef SingularDecomposition<MPStateComponent, MPStateComponent> AMatSVD;  // avoid typing...

#include "density.cc"

#endif

