// -*- C++ -*- $Id$

#if !defined(PRODUCTBASIS_H_JHCIUWHFIUH98743Y9843YP9)
#define PRODUCTBASIS_H_JHCIUWHFIUH98743Y9843YP9

#include "basis.h"
#include "tensor.h"
#include "linearalgebra/matrix.h"
//#include "linearalgebra/directproduct.h"
#include "pstream/pstream.h"

namespace Tensor
{

// This class maps the tensor product of LeftBasisType and RightBasisType onto a VectorBasis.
// LeftBasisType and RightBasisType must satisfy the VectorBasis interface.
template <typename BLType, typename BRType = BLType>
class ProductBasis;

template <typename BLType, typename BRType>
PStream::opstream& operator<<(PStream::opstream& out, ProductBasis<BLType, BRType> const& B);

template <typename BLType, typename BRType>
PStream::ipstream& operator>>(PStream::ipstream& in, ProductBasis<BLType, BRType>& B);

template <typename BLType, typename BRType>
class ProductBasis : public VectorBasis
{
   public:
      typedef BLType LeftBasisType;
      typedef BRType RightBasisType;
      typedef VectorBasis BasisType;

      // in the 'forward' direction, this maps a pair (s1, s2) onto a container of target subspaces
      typedef int target_type;
      // in the opposite direction, we map a subspace into a pair of original subspaces
      typedef std::pair<int, int> source_type;

   private:
      typedef std::list<target_type> TargetListType;
      typedef LinearAlgebra::Vector<source_type>    SourceListType;

   public:
      typedef TargetListType::const_iterator const_iterator;

      ProductBasis() {}

      ProductBasis(LeftBasisType Basis1_, RightBasisType Basis2_);

      ProductBasis(LeftBasisType Basis1_, RightBasisType Basis2_, QuantumNumber const& q);

      const_iterator begin(int s1, int s2) const { return TransformData(s1,s2).begin(); }
      const_iterator end(int s1, int s2) const { return TransformData(s1,s2).end(); }

      const_iterator begin(source_type const& s12) const { return this->begin(s12.first, s12.second); }
      const_iterator end(source_type const& s12) const { return this->end(s12.first, s12.second); }

      // the reverse mapping, from target subspace to source pair
      source_type rmap(int s) const { return ReverseMapping[s]; }

      // thunks the labels to the QuantumNumber types and returns the coefficient.
      double Coefficient(QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
			 int s1prime, int s2prime, int sprime, int s1, int s2, int s) const;
 
      LeftBasisType const& LeftBasis() const { return B1; }
      RightBasisType const& RightBasis() const { return B2; }

      int LeftSize() const { return B1.size(); }
      int RightSize() const { return B2.size(); }

      int LeftDimension(int s) const { return B1.Dimension(ReverseMapping[s].first); }
      int RightDimension(int s) const { return B2.Dimension(ReverseMapping[s].second); }

      VectorBasis const& Basis() const { return *this; }

   private:
      VectorBasis& Basis() { return *this; }

      LeftBasisType B1;
      RightBasisType B2;

      LinearAlgebra::Matrix<TargetListType> TransformData;
      SourceListType ReverseMapping;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ProductBasis const& B);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ProductBasis& B);
};

template <>
class ProductBasis<SimpleBasis, SimpleBasis> : public SimpleBasis
{
   public:
      typedef SimpleBasis LeftBasisType;
      typedef SimpleBasis RightBasisType;
      typedef SimpleBasis BasisType;

      // in the 'forward' direction, this maps a pair (s1, s2) onto a container of target subspaces
      typedef int target_type;
      // in the opposite direction, we map a subspace into a pair of original subspaces
      typedef std::pair<int, int> source_type;

   private:
      typedef std::list<target_type> TargetListType;
      typedef LinearAlgebra::Vector<source_type>    SourceListType;

   public:
      typedef TargetListType::const_iterator const_iterator;

      ProductBasis() {}

      ProductBasis(LeftBasisType Basis1_, RightBasisType Basis2_);

      // A product basis that projects onto a single quantum number only
      ProductBasis(LeftBasisType Basis1_, RightBasisType Basis2_, QuantumNumber const& q);

      // Named constructor for making a basis for a lower-triangular operator that projects onto some irreducible
      // component in the final element.  This is the same as a normal product basis between Basis1_ and Basis2_,
      // except that for the last element in Basis1_ and Basis2_, we take only the single projection q.
      static
      ProductBasis MakeTriangularProjected(LeftBasisType Basis1_, RightBasisType Basis2_, QuantumNumber const& q);

      const_iterator begin(int s1, int s2) const { return TransformData(s1,s2).begin(); }
      const_iterator end(int s1, int s2) const { return TransformData(s1,s2).end(); }

      const_iterator begin(source_type const& s12) const { return this->begin(s12.first, s12.second); }
      const_iterator end(source_type const& s12) const { return this->end(s12.first, s12.second); }

      // the reverse mapping, from target subspace to source pair
      source_type rmap(int s) const { return ReverseMapping[s]; }

      // thunks the labels to the QuantumNumber types and returns the coefficient.
      double Coefficient(QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
			 int s1prime, int s2prime, int sprime, int s1, int s2, int s) const;
 
      LeftBasisType const& LeftBasis() const { return B1; }
      RightBasisType const& RightBasis() const { return B2; }

      int LeftSize() const { return B1.size(); }
      int RightSize() const { return B2.size(); }

      SimpleBasis const& Basis() const { return *this; }

   private:
      ProductBasis(QuantumNumbers::SymmetryList const& sl) : SimpleBasis(sl) {}

      SimpleBasis& base() { return *this; }

      LeftBasisType B1;
      RightBasisType B2;

      LinearAlgebra::Matrix<TargetListType> TransformData;
      SourceListType ReverseMapping;

   friend PStream::opstream& operator<< <>(PStream::opstream& out, ProductBasis const& B);
   friend PStream::ipstream& operator>> <>(PStream::ipstream& in, ProductBasis& B);
};

template <typename BLType, typename BRType>
inline
ProductBasis<BLType, BRType>
MakeProductBasis(BLType const& BL, BRType const& BR)
{
   return ProductBasis<BLType, BRType>(BL, BR);
}

template <typename BLType, typename BRType>
inline
ProductBasis<BLType, BRType>
MakeProductBasis(BLType const& BL, BRType const& BR, QuantumNumber const& q)
{
   return ProductBasis<BLType, BRType>(BL, BR, q);
}

template <typename BLType, typename BRType>
inline double
ProductBasis<BLType, BRType>::Coefficient(QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
					  int s1prime, int s2prime, int sprime, int s1, int s2, int s) const
{
   return  tensor_coefficient(k1,             k2,             k,
			     B1.qn(s1prime), B2.qn(s2prime), this->qn(sprime),
			     B1.qn(s1),      B2.qn(s2),      this->qn(s));
}

template <typename BL1Type, typename BR1Type, typename BL2Type, typename BR2Type>
inline double 
tensor_coefficient(ProductBasis<BL1Type, BR1Type> const& P1, ProductBasis<BL2Type, BR2Type> const& P2,
		  QuantumNumber const& k1, QuantumNumber const& k2, QuantumNumber const& k,
		  int s1prime, int s2prime, int sprime, int s1, int s2, int s)
{
   //   TRACE(k1)(k2)(k)(P1.LeftBasis().qn(s1prime))( P1.RightBasis().qn(s2prime))(P1.qn(sprime));
   return  tensor_coefficient(k1,                         k2,                          k,
			     P1.LeftBasis().qn(s1prime), P1.RightBasis().qn(s2prime), P1.qn(sprime),
			     P2.LeftBasis().qn(s1),      P2.RightBasis().qn(s2),      P2.qn(s));
}
   
template <typename C1, typename BL1Type, typename BL2Type, 
	  typename C2, typename BR1Type, typename BR2Type,
	  typename ProductBasis1Type, typename ProductBasis2Type,
	  typename ProductFunctor>
IrredOperator<typename ProductFunctor::value_type, 
	      typename ProductBasis1Type::BasisType,
	      typename ProductBasis2Type::BasisType>
tensor_prod(IrredOperator<C1, BL1Type, BL2Type> const& ML, IrredOperator<C2, BR1Type, BR2Type> const& MR, 
	    QuantumNumber const& q, ProductBasis1Type const& PBasis1, ProductBasis2Type const& PBasis2,
	    ProductFunctor TensorProd)
{
   DEBUG_PRECONDITION_EQUAL(ML.Basis1(), PBasis1.LeftBasis());
   DEBUG_PRECONDITION_EQUAL(MR.Basis1(), PBasis1.RightBasis());
   DEBUG_PRECONDITION_EQUAL(ML.Basis2(), PBasis2.LeftBasis());
   DEBUG_PRECONDITION_EQUAL(MR.Basis2(), PBasis2.RightBasis());

   typedef IrredOperator<typename ProductFunctor::value_type, 
      typename ProductBasis1Type::BasisType,
      typename ProductBasis2Type::BasisType> ResultType;
   ResultType Result(PBasis1, PBasis2, q);

   typedef typename IrredOperator<C1, BL1Type, BL2Type>::MatrixType MLMatrixType;
   typedef typename MLMatrixType::row_vector_type MLRowType;
   typedef typename IrredOperator<C2, BR1Type, BR2Type>::MatrixType MRMatrixType;
   typedef typename MRMatrixType::row_vector_type MRRowType;
   typedef typename MLRowType::const_iterator MLIterType;
   typedef typename MRRowType::const_iterator MRIterType;

   for (int i1 = 0; i1 < ML.size1(); ++i1)
   {
      for (int i2 = 0; i2 < MR.size1(); ++i2)
      {
	 typename ProductBasis1Type::const_iterator TiIter, TiEnd = PBasis1.end(i1, i2);
	 for (TiIter = PBasis1.begin(i1, i2); TiIter != TiEnd; ++TiIter)
	 {
	    MLIterType j1End = ML.data().row_vector(i1).end();
	    for (MLIterType j1Iter = ML.data().row_vector(i1).begin(); j1Iter != j1End; ++j1Iter)
	    {
	       MRIterType j2End = MR.data().row_vector(i2).end();
	       for (MRIterType j2Iter = MR.data().row_vector(i2).begin(); j2Iter != j2End; ++j2Iter)
	       {
		  typename ProductBasis2Type::const_iterator TjIter, 
		     TjEnd = PBasis2.end(j1Iter.index(), j2Iter.index());
		  for (TjIter = PBasis2.begin(j1Iter.index(), j2Iter.index()); TjIter != TjEnd; ++TjIter)
		  {
		     double Coeff = tensor_coefficient(PBasis1, PBasis2, ML.TransformsAs(), MR.TransformsAs(), q,
						       i1, i2, *TiIter, j1Iter.index(), j2Iter.index(), *TjIter);

		     if (numerics::norm_2(Coeff) > 1E-10)
		     {
		        Result.data().assign_new(*TiIter, *TjIter, Coeff * TensorProd(j1Iter.data(), j2Iter.data()));
		     }
		  }
	       }
	    }
	 }
      }
   }
   return Result;
}

template <typename C1, typename BL1Type, typename BL2Type, 
	  typename C2, typename BR1Type, typename BR2Type,
	  typename ProductBasis1Type, typename ProductBasis2Type>
inline
IrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, 
	      typename ProductBasis1Type::BasisType,
	      typename ProductBasis2Type::BasisType>
tensor_prod(IrredOperator<C1, BL1Type, BL2Type> const& ML, IrredOperator<C2, BR1Type, BR2Type> const& MR, 
	    QuantumNumber const& q, ProductBasis1Type const& PBasis1, ProductBasis2Type const& PBasis2)
{
   return tensor_prod(ML, MR, q, PBasis1, PBasis2,  LinearAlgebra::DirectProduct<C1, C2>());
}

template <typename C1, typename BL1Type, typename BL2Type, 
	  typename C2, typename BR1Type, typename BR2Type,
	  typename ProductBasis1Type, typename ProductBasis2Type>
inline
IrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, 
	      typename ProductBasis<BL1Type, BL2Type>::BasisType,
	      typename ProductBasis<BR1Type, BR2Type>::BasisType>
tensor_prod(IrredOperator<C1, BL1Type, BL2Type> const& ML, IrredOperator<C2, BR1Type, BR2Type> const& MR, 
	    QuantumNumber const& q)
{
   return tensor_prod(ML, MR, q, ProductBasis<BL1Type, BL2Type>(ML.Basis1(), ML.Basis1()), 
		                 ProductBasis<BR1Type, BR2Type>(MR.Basis2(), MR.Basis2()));
}

//
// TensorProd
//
// A binary functor to calculate the tensor product of two operators, given
// the target quantum number and the product basis types as fixed.
//
// TL and TR are assumed to be IrredOperator's.
//

template <typename TL, typename TR, 
	  typename ProductBasis1Type = ProductBasis<typename TL::Basis1Type, typename TR::Basis1Type>,
	  typename ProductBasis2Type = ProductBasis<typename TL::Basis2Type, typename TR::Basis2Type> >
struct TensorProd
{
   TensorProd(QuantumNumber q_, ProductBasis1Type const& B1_, ProductBasis2Type const& B2_) 
      : q(q_), Basis1(B1_), Basis2(B2_) {}

   typedef typename TL::value_type CLType;
   typedef typename TR::value_type CRType;

   typedef IrredOperator<typename LinearAlgebra::DirectProduct<CLType, CRType>::value_type, 
			 typename ProductBasis1Type::BasisType,
			 typename ProductBasis2Type::BasisType>
      result_type;

   typedef result_type value_type;

   result_type operator()(TL const& A, TR const& B) const
   {
      return tensor_prod(A, B, q, Basis1, Basis2);
   }

   QuantumNumber q;
   ProductBasis1Type Basis1;
   ProductBasis2Type Basis2;
};

#if 0

template <typename C1, typename C2>
inline
ComplexIrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, VectorBasis>
tensor_prod(IrredOperator<C1, VectorBasis> const& M1, ComplexIrredOperator<C2, VectorBasis> const& M2, 
	    QuantumNumber const& q, ProductBasis const& PBasis)
{
   typedef ComplexIrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, VectorBasis> res_type;
   res_type Result(PBasis, q);
   Result.real() = tensor_prod(M1, M2.real(), q, PBasis);
   Result.imag() = tensor_prod(M1, M2.imag(), q, PBasis);
   return Result;
}

template <typename C1, typename C2>
inline
ComplexIrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, VectorBasis>
tensor_prod(ComplexIrredOperator<C1, VectorBasis> const& M1, IrredOperator<C2, VectorBasis> const& M2, 
	    QuantumNumber const& q, ProductBasis const& PBasis)
{
   typedef ComplexIrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, VectorBasis> res_type;
   res_type Result(PBasis, q);
   Result.real() = tensor_prod(M1.real(), M2, q, PBasis);
   Result.imag() = tensor_prod(M1.imag(), M2, q, PBasis);
   return Result;
}

template <typename C1, typename C2>
inline
ComplexIrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, VectorBasis>
tensor_prod(ComplexIrredOperator<C1, VectorBasis> const& M1, ComplexIrredOperator<C2, VectorBasis> const& M2, 
	    QuantumNumber const& q, ProductBasis const& PBasis)
{
   typedef ComplexIrredOperator<typename LinearAlgebra::DirectProduct<C1, C2>::value_type, VectorBasis> res_type;
   res_type Result(PBasis, q);
   Result.real() = tensor_prod(M1.real(), M2.real(), q, PBasis);
   Result.real() -= tensor_prod(M1.imag(), M2.imag(), q, PBasis);
   Result.imag() = tensor_prod(M1.real(), M2.imag(), q, PBasis);
   Result.imag() += tensor_prod(M1.imag(), M2.real(), q, PBasis);
   return Result;
}

inline
Operator tensor_prod(Operator const& M1, Operator const& M2, QuantumNumber const& q, ProductBasis const& PBasis)
{
   PRECONDITION(M1.GetSymmetryList() == M2.GetSymmetryList());

   Operator Result(PBasis, q);
   Result += tensor_prod(M1.Dense(), M2.Dense(), q, PBasis);
   Result += tensor_prod(M1.Dense(), M2.Sparse(), q, PBasis);
   Result += tensor_prod(M1.Sparse(), M2.Dense(), q, PBasis);
   Result += tensor_prod(M1.Sparse(), M2.Sparse(), q, PBasis);

   return Result;
}

#endif

} // namespace Tensor

#include "productbasis.cc"

#endif
