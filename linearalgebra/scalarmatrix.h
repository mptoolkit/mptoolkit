// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/scalarmatrix.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
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

/*
  scalarmatrix.h

  A diagonal matrix where the diagonal is a fixed number.  Ie, this
  represents a constant times the identity matrix.

  Created 2005-03-08 Ian McCulloch
*/

#if !defined(SCALARMATRIX_H_CJDHCJKHUIYT389Y987Y34897GYO8)
#define SCALARMATRIX_H_CJDHCJKHUIYT389Y987Y34897GYO8

#include "fixedvector.h"
#include "matrixiterators.h"
#include "matrixdirectproduct.h"
#include "matrixaddition.h"
#include "crtp_matrix.h"
#include "matrixmatrixmultiplication.h"

namespace LinearAlgebra
{

// tagged constructor for direct initialization
struct cdirect {};

template <typename T>
class ScalarMatrix : public MatrixBase<ScalarMatrix<T> >
{
   public:
      typedef T data_type;

      typedef is_mutable_proxy<data_type> proxy;
      typedef is_const_proxy<data_type> const_proxy;
      typedef is_immediate<data_type> immediate;

      typedef typename make_value<T>::type           value_type;
      typedef typename make_reference<T>::type       reference;
      typedef typename make_const_reference<T>::type const_reference;

      ScalarMatrix() : size_(0) {}

      explicit ScalarMatrix(size_type s1, const_reference value = T())
         : size_(s1), value_(value) { }

      ScalarMatrix(size_type s1, size_type s2, const_reference value = T())
         : size_(s1), value_(value) { DEBUG_CHECK_EQUAL(s1,s2); }

      explicit ScalarMatrix(size_type s1, reference value)
         : size_(s1), value_(value) { }

      ScalarMatrix(size_type s1, size_type s2, reference value)
         : size_(s1), value_(value) { DEBUG_CHECK_EQUAL(s1,s2); }

      template <typename U>
      ScalarMatrix(size_type size, U const& value, cdirect)
         : size_(size), value_(value) {}

      template <typename U>
      ScalarMatrix(size_type size, U& value, cdirect)
         : size_(size), value_(value) {}

      size_type size1() const { return size_; }
      size_type size2() const { return size_; }

      reference value() { return value_; }
      const_reference value() const { return value_; }

   private:
      size_type size_;
      data_type value_;
};

// interface

template <typename T>
struct interface<ScalarMatrix<T> >
{
   using value_type = typename make_value<T>::type;
   using type = Concepts::ScalarMatrix<value_type, ScalarMatrix<T>>;
};

// iterators

template <typename T>
struct Iterate<ScalarMatrix<T> >
{
   typedef ConstantIterator<T> iviter;
   typedef MatrixIterDiagonal<iviter> imiter;
   typedef MatrixDummyOuterIterator<imiter> result_type;
   typedef ScalarMatrix<T> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(imiter(iviter(x.size1(), x.value())));
   }
};

// GetMatrixElement

// we go to some lengths here to get this to work with value_type's that
// don't have a zero value.

template <typename T, typename Enable = void>
struct GetMatrixElement_ScalarMatrix {};

template <typename T>
struct GetMatrixElement_ScalarMatrix<
   T
 , typename boost::enable_if<has_zero<typename make_value<T>::type> >::type
>
{
   typedef ScalarMatrix<T> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef typename ScalarMatrix<T>::const_reference result_type;

   result_type operator()(first_argument_type x, size_type i, size_type j) const
   {
      DEBUG_PRECONDITION_COMPARE(i, <, x.size1());
      DEBUG_PRECONDITION_COMPARE(j, <, x.size2());
      if (i == j) return x.value();
      else return StaticZero<typename make_value<T>::type>::value;
   }
};

template <typename T>
struct GetMatrixElement_ScalarMatrix<
   T
 , typename boost::disable_if<has_zero<typename make_value<T>::type> >::type
>
{
   typedef ScalarMatrix<T> const& first_argument_type;
   typedef size_type second_argument_type;
   typedef size_type third_argument_type;
   typedef typename make_value_with_zero<T>::type result_type;

   result_type operator()(first_argument_type x, size_type i, size_type j) const
   {
      DEBUG_PRECONDITION_COMPARE(i, <, x.size1());
      DEBUG_PRECONDITION_COMPARE(j, <, x.size2());
      return (i == j) ? x.value() : zero<result_type>();
   }
};

template <typename T>
struct GetMatrixElement<ScalarMatrix<T> > : GetMatrixElement_ScalarMatrix<T> {};

// disable the mutable version of GetMatrixElement
template <typename T>
struct GetMatrixElement<ScalarMatrix<T>&> {};

//
// expressions
//

// Conj

template <typename T>
struct Conj<ScalarMatrix<T> >
{
   typedef ScalarMatrix<T> const& argument_type;
   typedef ScalarMatrix<typename result_value<Conj<T> >::type> result_type;
   typedef typename is_identity<Conj<T> >::type identity;

   result_type operator()(argument_type x) const
   {
      return result_type(conj(x.value()));
   }
};

template <typename T>
struct Herm<ScalarMatrix<T> >
{
   typedef ScalarMatrix<T> const& argument_type;
   typedef ScalarMatrix<typename result_value<Conj<T> >::type> result_type;
   typedef typename is_identity<Conj<T> >::type identity;

   result_type operator()(argument_type x) const
   {
      return result_type(x.size1(), conj(x.value()));
   }
};

// Multiply

#if 0
template <typename T, typename U>
struct Multiply<ScalarMatrix<T>&, U, typename boost::enable_if<is_defined<Multiply<T&, U>>>::type>
{
   typedef ScalarMatrix<T>& result_type;
   typedef ScalarMatrix<T>& first_argument_type;
   typedef U second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      multiply(x.value(), y);
      return x;
   }
};
#endif

template <typename T, typename RHS, typename S1, typename U1, typename RHSi>
struct MultiplyInterface<ScalarMatrix<T>&, RHS, Concepts::ScalarMatrix<S1, U1>, AnyScalar<RHSi>>
{
   typedef ScalarMatrix<T>& result_type;
   typedef ScalarMatrix<T>& first_argument_type;
   typedef RHS const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      DEBUG_CHECK_EQUAL(size2(x), size1(y));
      multiply(x.value(), y);
      return x;
   }
};

// Addition

template <typename LHS, typename RHS,
          typename S1, typename U1,
          typename S2, typename U2>
struct AddInterface<LHS&, RHS, Concepts::ScalarMatrix<S1, U1>, Concepts::ScalarMatrix<S2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_CHECK_EQUAL(size2(x), size1(y));
      x.value() += y.value();
   }
};

template <typename LHS, typename RHS,
          typename S1, typename U1,
          typename S2, typename U2>
struct SubtractInterface<LHS&, RHS, Concepts::ScalarMatrix<S1, U1>, Concepts::ScalarMatrix<S2, U2>>
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef RHS const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_CHECK_EQUAL(size1(x), size1(y));
      x.value() -= y.value();
   }
};

// Transform

template <typename T, typename F>
struct TransformMatrix<ScalarMatrix<T>, F>
{
   typedef ScalarMatrix<T> const& argument_type;
   typedef argument_type first_argument_type;
   typedef F const& second_argument_type;

   typedef ScalarMatrix<typename make_value<typename F::result_type>::type> result_type;

   result_type operator()(argument_type e) const
   {
      return result_type(e.size1(), F()(e.value()), cdirect());
   }

   result_type operator()(argument_type e, second_argument_type f) const
   {
      return result_type(e.size1(), f(e.value()), cdirect());
   }
};

template <typename T, typename F>
struct TransformMatrix<ScalarMatrix<T>&, F>
{
   typedef ScalarMatrix<T>& argument_type;
   typedef argument_type first_argument_type;
   typedef F const& second_argument_type;

   typedef ScalarMatrix<typename F::result_type> result_type;

   result_type operator()(argument_type e) const
   {
      return result_type(e.size1(), F()(e.value()), cdirect());
   }

   result_type operator()(argument_type e, second_argument_type f) const
   {
      return result_type(e.size1(), f(e.value()), cdirect());
   }
};

// BinaryTransform

template <typename S, typename T, typename F, typename Sv, typename Tv, typename Si, typename Ti>
struct BinaryTransformMatrixSemiregular<ScalarMatrix<S>, ScalarMatrix<T>, F,
                                        Concepts::ScalarMatrix<Sv, Si>,
                                        Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef ScalarMatrix<S> const& first_argument_type;
   typedef ScalarMatrix<T> const& second_argument_type;
   typedef F const& third_argument_type;

   typedef ScalarMatrix<typename make_value<typename F::result_type>::type> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return result_type(x.size1(), F()(x.value(), y.value()), cdirect());
   }

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type f) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return result_type(x.size1(), f(x.value(), y.value()), cdirect());
   }
};

template <typename S, typename T, typename Sv, typename Tv, typename Si, typename Ti>
struct BinaryTransformMatrixSemiregular<ScalarMatrix<S>, ScalarMatrix<T>, Addition<Sv, Tv>,
                                        Concepts::ScalarMatrix<Sv, Si>,
                                        Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef Addition<Sv, Tv> F;
   typedef ScalarMatrix<S> const& first_argument_type;
   typedef ScalarMatrix<T> const& second_argument_type;
   typedef F const& third_argument_type;

   typedef ScalarMatrix<typename make_value<typename F::result_type>::type> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return result_type(x.size1(), F()(x.value(), y.value()), cdirect());
   }

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type f) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return result_type(x.size1(), f(x.value(), y.value()), cdirect());
   }
};

template <typename S, typename T, typename Sv, typename Tv, typename Si, typename Ti>
struct BinaryTransformMatrixSemiregular<ScalarMatrix<S>, ScalarMatrix<T>, Subtraction<Sv, Tv>,
                                        Concepts::ScalarMatrix<Sv, Si>,
                                        Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef Subtraction<Sv, Tv> F;
   typedef ScalarMatrix<S> const& first_argument_type;
   typedef ScalarMatrix<T> const& second_argument_type;
   typedef F const& third_argument_type;

   typedef ScalarMatrix<typename make_value<typename F::result_type>::type> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return result_type(x.size1(), F()(x.value(), y.value()), cdirect());
   }

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type f) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return result_type(x.size1(), f(x.value(), y.value()), cdirect());
   }
};

// multiplication

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<S, T, F,
                                  Concepts::ScalarMatrix<Sv, Si>,
                                  Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef typename is_commutative<F>::type commutative;

   typedef typename F::result_type FwdResult;
   typedef typename make_value<FwdResult>::type result_value_type;

   typedef ScalarMatrix<result_value_type> result_type;
   using first_argument_type = S;
   using second_argument_type = T;

   result_type operator()(S const& x, T const& y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return result_type(x.size1(), F()(x.value(), y.value()), cdirect());
   }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<S, T, F,
                                  Concepts::ScalarMatrix<Sv, Si>,
                                  Concepts::AnyMatrix<Tv,Ti>>
{
   typedef typename is_commutative<F>::type commutative;

   typedef ScalarMatrixMultiplication<Sv, T, F> Fwd;
   typedef typename Fwd::result_type result_type;
   typedef ScalarMatrix<S> const& first_argument_type;
   typedef T const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return Fwd()(x.value(), y);
   }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<S, T, F,
                                  Concepts::AnyMatrix<Sv, Si>,
                                  Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef MatrixScalarMultiplication<S, Tv, F> Fwd;
   typedef typename Fwd::result_type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return Fwd()(x, y.value());
   }
};

//
// DirectProduct
//

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixDirectProduct<ScalarMatrix<S>, ScalarMatrix<T>, F,
                           Concepts::ScalarMatrix<Sv, Si>,
                           Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef typename F::result_type ValType;
   typedef typename make_value<ValType>::type result_value_type;
   typedef ScalarMatrix<result_value_type> result_type;

   typedef ScalarMatrix<S> const& first_argument_type;
   typedef ScalarMatrix<T> const& second_argument_type;
   typedef F const& third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.size1() * y.size1(), F()(x.value(), y.value()), cdirect());
   }

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type f) const
   {
      return result_type(x.size1() * y.size1(), f(x.value(), y.value()), cdirect());
   }
};

// InnerProd

template <typename S, typename T>
struct InnerProd<ScalarMatrix<S>, ScalarMatrix<T>>
{
   typedef ScalarMatrix<S> const& first_argument_type;
   typedef ScalarMatrix<T> const& second_argument_type;
   typedef typename result_value<InnerProd<S, T> >::type result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return inner_prod(x.value(), y.value()) * result_type(x.size1());
   }
};

template <typename S, typename T,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<S, ScalarMatrix<T>, InnerProd<Sv, Tv>,
                       Concepts::AnyMatrix<Sv, Si>,
                       Concepts::ScalarMatrix<Tv, Ti>>
{
   typedef S const& first_argument_type;
   typedef ScalarMatrix<T> const& second_argument_type;
   typedef typename result_value<InnerProd<Sv, Tv> >::type result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return inner_prod(trace(x), y.value());
   }
};

template <typename S, typename T,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<ScalarMatrix<S>, T, InnerProd<Sv, Tv>,
                       Concepts::ScalarMatrix<Sv, Si>,
                       Concepts::AnyMatrix<Tv, Ti>>
{
   typedef ScalarMatrix<S> const& first_argument_type;
   typedef T const& second_argument_type;
   typedef typename result_value<InnerProd<Sv, Tv> >::type result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return inner_prod(x.value(), trace(y));
   }
};



} // namespace LinearAlgebra

#endif
