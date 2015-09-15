/* -*- C++ -*- $Id$

  diagonalmatrix.h

  Created 2009-03-17 Ian McCulloch
*/

#if !defined(DIAGONALMATRIX_H_CJDHCJKHUIYT389Y987Y34897GYO8)
#define DIAGONALMATRIX_H_CJDHCJKHUIYT389Y987Y34897GYO8

#include "fixedvector.h"
#include "matrixiterators.h"
#include "matrixdirectproduct.h"
#include "matrixaddition.h"
#include "crtp_matrix.h"
#include "matrixmatrixmultiplication.h"

namespace LinearAlgebra
{

template <typename T, typename DiagonalType = Vector<T> >
class DiagonalMatrix : public MatrixBase<DiagonalMatrix<T> >
{
   public:
      typedef T data_type;

      typedef is_mutable_proxy<data_type> proxy;
      typedef is_const_proxy<data_type> const_proxy;
      typedef is_immediate<data_type> immediate; 

      typedef typename make_value<T>::type           value_type;
      typedef typename make_reference<T>::type       reference;
      typedef typename make_const_reference<T>::type const_reference;
      typedef DiagonalType diagonal_type;

      DiagonalMatrix() {}

      template <typename U>
      DiagonalMatrix(DiagonalMatrix<U> const& Other)
	 : data_(Other.diagonal()) {}

      DiagonalMatrix(size_type s) : data_(s) { }

      DiagonalMatrix(size_type s1, size_type s2, const_reference value = T())
         : data_(s1, value) { DEBUG_CHECK_EQUAL(s1,s2); }

      DiagonalMatrix(size_type s1, size_type s2, reference value)
         : data_(s1, value) { DEBUG_CHECK_EQUAL(s1,s2); }

      // initialization from the diagonal vector
      explicit DiagonalMatrix(diagonal_type const& diag)
	 : data_(diag) {}
      
      size_type size1() const { return size(data_); }
      size_type size2() const { return size(data_); }

      DiagonalMatrix& operator=(DiagonalMatrix const& x)
      {
	 data_ = x.data_;
	 return *this;
      }

      diagonal_type const& diagonal() const { return data_; }
      diagonal_type& diagonal() { return data_; }

   private:
      diagonal_type data_;
};

// interface

template <typename T>
struct interface<DiagonalMatrix<T> >
{
   typedef typename make_value<T>::type value_type;
   typedef DIAGONAL_MATRIX(value_type, DiagonalMatrix<T>) type;
};

// iterators

template <typename T>
struct Iterate<DiagonalMatrix<T> >
{
   typedef typename const_iterator<typename DiagonalMatrix<T>::diagonal_type>::type diag_iter;
   typedef MatrixIterDiagonal<diag_iter> imiter;
   typedef MatrixDummyOuterIterator<imiter> result_type;
   typedef DiagonalMatrix<T> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(imiter(iterate(x.diagonal())));
   }
};

template <typename T>
struct Iterate<DiagonalMatrix<T>&>
{
   typedef typename iterator<typename DiagonalMatrix<T>::diagonal_type>::type diag_iter;
   typedef MatrixIterDiagonal<diag_iter> imiter;
   typedef MatrixDummyOuterIterator<imiter> result_type;
   typedef DiagonalMatrix<T>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(imiter(iviter(iterate(x.diagonal()))));
   }
};

template <typename T>
struct IterateDiagonal {};

template <typename T>
struct IterateDiagonal<DiagonalMatrix<T> >
{
   typedef typename DiagonalMatrix<T>::diagonal_type diagonal_type;
   typedef typename const_iterator<diagonal_type>::type result_type;
   typedef DiagonalMatrix<T> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return iterate(x.diagonal());
   }
};

template <typename T>
struct IterateDiagonal<DiagonalMatrix<T>&>
{
   typedef typename DiagonalMatrix<T>::diagonal_type diagonal_type;
   typedef typename iterator<diagonal_type>::type result_type;
   typedef DiagonalMatrix<T>& argument_type;
   result_type operator()(argument_type x) const
   {
      return iterate(x.diagonal());
   }
};

// GetMatrixElement

// we go to some lengths here to get this to work with value_type's that
// don't have a zero value.

template <typename T, typename Enable = void>
struct GetMatrixElement_DiagonalMatrix {};

template <typename T>
struct GetMatrixElement_DiagonalMatrix<
   T
 , typename boost::enable_if<has_zero<typename make_value<T>::type> >::type
>
{
   typedef DiagonalMatrix<T> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef typename DiagonalMatrix<T>::const_reference result_type;

   result_type operator()(first_argument_type x, size_type i, size_type j) const
   {
      DEBUG_PRECONDITION_COMPARE(i, <, x.size1());
      DEBUG_PRECONDITION_COMPARE(j, <, x.size2());
      if (i == j) return x.diagonal()[i];
      else return StaticZero<typename make_value<T>::type>::value;
   }
};

template <typename T>
struct GetMatrixElement_DiagonalMatrix<
   T
 , typename boost::disable_if<has_zero<typename make_value<T>::type> >::type
>
{
   typedef DiagonalMatrix<T> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef typename make_value_with_zero<T>::type result_type;

   result_type operator()(first_argument_type x, size_type i, size_type j) const
   {
      DEBUG_PRECONDITION_COMPARE(i, <, x.size1());
      DEBUG_PRECONDITION_COMPARE(j, <, x.size2());
      return (i == j) ? x.diagonal()[i] : zero<result_type>();
   }
};

template <typename T>
struct GetMatrixElement<DiagonalMatrix<T> > : GetMatrixElement_DiagonalMatrix<T> {};

// mutable version of GetMatrixElement
template <typename T>
struct GetMatrixElement<DiagonalMatrix<T>&>
{
   typedef DiagonalMatrix<T>& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef typename DiagonalMatrix<T>::reference result_type;

   result_type operator()(first_argument_type x, size_type i, size_type j) const
   {
      DEBUG_PRECONDITION_COMPARE(i, <, x.size1());
      DEBUG_PRECONDITION_COMPARE(j, <, x.size2());
      DEBUG_PRECONDITION_EQUAL(i,j);
      return x.diagonal()[i];
   }
};

//
// expressions
//

// Multiply

template <typename T, typename RHS, typename S1, typename U1, typename RHSi>
struct MultiplyInterface<DiagonalMatrix<T>&, RHS, DIAGONAL_MATRIX(S1, U1), AnyScalar<RHSi> >
{
   typedef DiagonalMatrix<T>& result_type;
   typedef DiagonalMatrix<T>& first_argument_type;
   typedef RHS const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      multiply(x.diagonal(), y);
      return x;
   }
};

// multiplication

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<DiagonalMatrix<S>, DiagonalMatrix<T>, F,
                                  DIAGONAL_MATRIX(Sv, Si), DIAGONAL_MATRIX(Tv, Ti)>
{
   typedef typename is_commutative<F>::type commutative;

   typedef typename F::result_type FwdResult;
   typedef typename make_value<FwdResult>::type result_value_type;

   typedef DiagonalMatrix<result_value_type> result_type;
   typedef DiagonalMatrix<S> const& first_argument_type;
   typedef DiagonalMatrix<T> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return result_type(transform(x.diagonal(), y.diagonal(), F()));
   }
};

//
// DirectProduct
//

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixDirectProduct<DiagonalMatrix<S>, DiagonalMatrix<T>, F, 
                           DIAGONAL_MATRIX(Sv, Si), DIAGONAL_MATRIX(Tv, Ti)>
{
   typedef typename F::result_type ValType;
   typedef typename make_value<ValType>::type result_value_type;
   typedef DiagonalMatrix<result_value_type> result_type;

   typedef DiagonalMatrix<S> const& first_argument_type;
   typedef DiagonalMatrix<T> const& second_argument_type;
   typedef F const& third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(direct_product(x.diagonal(), y.diagonal(), F()));
   }

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type f) const
   {
      return result_type(direct_product(x.diagonal(), y.diagonal(), f));
   }
};

// InnerProd

template <typename S, typename T>
struct InnerProd<DiagonalMatrix<S>, DiagonalMatrix<T> >
{
   typedef DiagonalMatrix<S> const& first_argument_type;
   typedef DiagonalMatrix<T> const& second_argument_type;
   typedef typename result_value<InnerProd<S, T> >::type result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return inner_prod(x.diagonal(), y.diagonal());
   }
};

} // namespace LinearAlgebra

#include "matrixproductdiagonal.h"

#endif
