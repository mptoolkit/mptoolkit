/* -*- C++ -*- $Id$

  Adds a VectorBasis to the IrredTensor.
*/

#if !defined(BLOCKOPERATOR_H_DCKJHFAUOIFHIU4YU58943Y9PFWE)
#define BLOCKOPERATOR_H_DCKJHFAUOIFHIU4YU58943Y9PFWE

#include "tensor/tensor.h"
#include "tensor/tensorsum.h"
#include "tensor/tensorproduct.h"
#include "linearalgebra/scalarmatrix.h"

using namespace Tensor;

typedef IrredTensor<LinearAlgebra::Matrix<std::complex> > DenseOperator;

template <typename T, typename S = typename IrredTensor<T>::MatrixType>
class BlockOperator
{
   public:
      typedef IrredTensor<T, S> base_type;
      typedef typename base_type::value_type value_type;

      typedef VectorBasis basis1_type;
      typedef VectorBasis basis2_type;

      BlockOperator();

      BlockOperator(VectorBasis const& b1, VectorBasis const& b2, QuantumNumber const& q);

      BlockOperator(VectorBasis const& b1, VectorBasis const& b2, base_type const& b);

      VectorBasis const& Basis1() const { return Basis1_; }
      VectorBasis const& Basis2() const { return Basis2_; }

      QuantumNumber const& TransformsAs() const { return SparsePart_.TransformsAs(); }

      BlockOperator& operator+=(BlockOperator const& x);
      BlockOperator& operator-=(BlockOperator const& x);

      BlockOperator& operator*=(double a);

      SymmetryList GetSymmetryList() const { return Base_.GetSymmetryList(); }

      size_type size1() const { return Basis1_.size(); }
      size_type size2() const { return Basis2_.size(); }

      QuantumNumber const& qn1(size_type i) const { return Basis1_[i]; }
      QuantumNumber const& qn2(size_type j) const { return Basis2_[j]; }

      base_type& base() { return Base_; }
      base_type const& base() const { return Base_; }

      operator base_type&() { return Base_; }
      operator base_type const&() const { return Base_; }

   private:
      VectorBasis Basis1_, Basis2_;
      base_type Base_;
};

namespace LinearAlgebra
{

template <typename T>
struct interface<BlockOperator<T> >
{
   typedef void type;
};

// comparison

template <typename T, typename S>
struct Equal<BlockOperator<T, S>, BlockOperator<T, S> >
{
   typedef BlockOperator<T, S> const& first_argument_type;
   typedef BlockOperator<T, S> const& second_argument_type;
   typedef bool result_type;

   Equal(double tol =  LinearAlgebra::default_tolerance()) : tol_(tol) {}

   bool operator()(first_argument_type x, second_argument_type y) const
   {
      return x.TransformsAs() == y.TransformsAs() && norm_frob(x-y) <= tol_;
   }

   double tol_;
};

// arithmetic

template <typename T, typename S>
struct Negate<BlockOperator<T, S> >
{
   typedef BlockOperator<T, S> const& argument_type;
   typedef BlockOperator<T< S> result_type;

   bool operator()(argument_type x) const
   {
      return BlockOperator(x.Basis1(), x.Basis2(), -x.base());
   }
};

template <typename T, typename S>
struct Addition<BlockOperator<T, S>, BlockOperator<T, S> >
{
   typedef BlockOperator<T, S> const& first_argument_type;
   typedef BlockOperator<T, S> const& second_argument_type;
   typedef BlockOperator<T, S> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      return result_type(x.Basis1(), x.Basis2(), x.base()+y.base());
   }
};

template <>
struct Subtraction<BlockOperator, BlockOperator>
{
   typedef BlockOperator const& first_argument_type;
   typedef BlockOperator const& second_argument_type;
   typedef BlockOperator result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.Basis1(), y.Basis1());
      DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.Basis2());
      return result_type(x.Basis1(), x.Basis2(), x.base()-y.base());
   }
};

template <typename U, typename T, typename S>
struct Multiplication<U, BlockOperator<T, S> >
{
   typedef U first_argument_type;
   typedef BlockOperator<T, S> const& second_argument_type;
   typedef BlockOperator<T< S> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return BlockOperator(y.Basis1(), y.Basis2(), x*y.base());
   }
};

template <typename T, typename S, typename U>
struct Multiplication<BlockOperator, U>
{
   typedef BlockOperator<T, S> const& first_argument_type;
   typedef U second_argument_type;
   typedef BlockOperator<T, S> result_type;

   result_type operator()(first_argument_type y, second_argument_type x) const
   {
      return BlockOperator(y.Basis1(), y.Basis2(), y.base()*x);
   }
};

// complex conjugation

template <typename T, typename S>
struct Conj<BlockOperator<T, S> >
{
   typedef BlockOperator<T, S> const& argument_type;
   typedef BlockOperator<T, S> result_type;

   bool operator()(argument_type x) const
   {
      return BlockOperator(x.Basis1(), x.Basis2(), conj(x.base()));
   }
};
   
// Hermitian conjugation

template <typename T, typename S>
struct Herm<BlockOperator>
{
   typedef BlockOperator<T, S> const& argument_type;
   typedef HermitianProxy<BlockOperator<T, S> > result_type;

   bool operator()(argument_type x) const
   {
      return result_type(x);
   }
};

// Trace

template <typename T, typename S>
struct Trace<BlockOperator<T, S> > : Trace<IrredTensor<T, S> > {};

// norm_frob

template <typename T, typename S>
struct NormFrobSq<BlockOperator<T, S> > : NormFrobSq<IrredTensor<T, S> > {};

// inner_prod

template <typename T1, typename S1, typename T2, typename S2>
struct InnerProd<BlockOperator<T1, S1>, BlockOperator<T2, S2> > 
   : InnerProd<IrredTensor<T1, S1>, IrredTensor<T2, S2> > {};


// scalar_prod

template <typename T1, typename T2, typename S, typename Functor>
struct ScalarProd<BlockOperator<T1, S>, BlockOperator<T2, S>, Functor>
{
   typedef BlockOperator<T1, S> const& first_argument_type;
   typedef BlockOperator<T2, S> const& second_argument_type;
   typedef BlockOperator<typename result_value<Functor>::type> result_type;

   ScalarProd(Functor f = Functor()) : f_(f) {}
   
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.Basis1(), y.Basis2(), scalar_prod(x.base(), y.base(), f_));
   }

   Functor f_;
};

template <typename T1, typename T2, typename S, typename Functor>
struct ScalarProd<BlockOperator<T1, S>, HermitianProxy<BlockOperator<T2, S> >, Functor>
{
   typedef BlockOperator<T1, S> const& first_argument_type;
   typedef HermitianProxy<BlockOperator<T2, S> > const& second_argument_type;
   typedef BlockOperator<typename result_value<Functor>::type> result_type;

   ScalarProd(Functor f = Functor()) : f_(f) {}
   
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.Basis1(), y.Basis2(), scalar_prod(x.base(), herm(y.base().base()), f_));
   }

   Functor f_;
};

template <typename T1, typename T2, typename S, typename Functor>
struct ScalarProd<HermitianProxy<BlockOperator<T1, S> >, BlockOperator<T2, S>, Functor>
{
   typedef HermitianProxy<BlockOperator<T1, S> > const& first_argument_type;
   typedef BlockOperator<T2, S> const& second_argument_type;
   typedef BlockOperator<typename result_value<Functor>::type> result_type;

   ScalarProd(Functor f = Functor()) : f_(f) {}
   
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.Basis1(), y.Basis2(), scalar_prod(herm(x.base().base()), y.base(), f_));
   }

   Functor f_;
};

// prod

template <typename T1, typename S1, typename T2, typename S2, typename Nest>
struct IrredProd<BlockOperator<T1, S1>, BlockOperator<T2, S2>, Nest>
{
   typedef BlockOperator<T1, S1> const& first_argument_type;
   typedef BlockOperator<T2, S2> const& second_argument_type;
   typedef Nest third_argument_type;
   typedef BlockOperator<typename result_value<Nest>::type> result_type;

   result_type operator()(first_argument_type x, second_argument_type y,
                          Nest F) const
   {
      return result_type(x.Basis1(), y.Basis2(), prod(x.base(), y.base(), F));
   }

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.Basis1(), y.Basis2(), prod(x.base(), y.base()));
   }
};

template <typename T1, typename S1, typename T2, typename S2>
struct Multiplication<BlockOperator<T1, S1>, BlockOperator<T2, S2> >
   : IrredProd<BlockOperator<T1, S1>, BlockOperator<T2, S2> > {};

} // namespace LinearAlgebra

// triple_prod

template <typename T1, typename S1, typename T2, typename S2, typename T3, typename S3>
BlockOperator triple_prod(BlockOperator<T1, S1> const& x, 
			  BlockOperator<T2, S2> const& E,
			  HermitianProxy<BlockOperator<T3, S3> > y);

template <typename T1, typename S1, typename T2, typename S2, typename T3, typename S3>
BlockOperator triple_prod(BlockOperator<T1, S1> const& x, 
			  BlockOperator<T2, S2> const& E,
			  HermitianProxy<BlockOperator const> y,
			  QuantumNumber qxy,
			  QuantumNumber qEp);

template <typename T1, typename S1, typename T2, typename S2, typename T3, typename S3>
BlockOperator triple_prod(HermitianProxy<BlockOperator const> x, 
			  BlockOperator const& E,
			  BlockOperator const& y);

template <typename T1, typename S1, typename T2, typename S2, typename T3, typename S3>
BlockOperator triple_prod(HermitianProxy<BlockOperator const> x, 
			  BlockOperator const& E,
			  BlockOperator const& y,
			  QuantumNumber qxy,
			  QuantumNumber qEp);

BlockOperator adjoint(BlockOperator const& x);
BlockOperator inv_adjoint(BlockOperator const& x);

BlockOperator tensor_prod(BlockOperator const& x, BlockOperator const& y,
			  VectorProductBasis const& b1,
			  VectorProductBasis const& b2,
			  QuantumNumber const& q = QuantumNumber());

BlockOperator tensor_sum(BlockOperator const& x, BlockOperator const& y,
			 VectorSumBasis const& b1, VectorSumBasis const& b2);

BlockOperator tensor_row_sum(BlockOperator const& x, BlockOperator const& y,
			     VectorSumBasis const& b1);

BlockOperator tensor_col_sum(BlockOperator const& x, BlockOperator const& y,
			     VectorSumBasis const& b2);

#if 0
template <typename FwdIter>
typename std::iterator_traits<FwdIter>::value_type
tensor_accumulate(FwdIter first, FwdIter last,
		  VectorSumBasis const& B1, 
		  VectorSumBasis const& B2);

template <typename FwdIter>
typename std::iterator_traits<FwdIter>::value_type
tensor_row_accumulate(FwdIter first, FwdIter last, VectorSumBasis const& B2);

template <typename FwdIter>
typename std::iterator_traits<FwdIter>::value_type
tensor_col_accumulate(FwdIter first, FwdIter last, VectorSumBasis const& B1);
#endif

#include "blockoperator.cc"

#endif
