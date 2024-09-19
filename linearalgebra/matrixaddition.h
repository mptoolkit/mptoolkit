// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixaddition.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
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

/*
  matrixaddition.h

  A proxy class & iterators for binary expressions of matrices

  Created 2005-03-15 Ian McCulloch
*/

#if !defined(MPTOOLKIT_LINEARALGEBRA_MATRIXADDITION_H)
#define MPTOOLKIT_LINEARALGEBRA_MATRIXADDITION_H

#include "interface.h"
#include "matrixoperationsbase.h"
#include "matrixbinaryiterator.h"
#include <boost/mpl/assert.hpp>

namespace LinearAlgebra
{

template <typename T1, typename T2>
class MatrixAdditionProxy;

// It isn't straightforward to determine the interface concept of
// a MatrixAdditionProxy, so we do a bit of metaprogramming.
// If the matrices are dense, then the interface is a DenseMatrix.
// Otherwise we just take the interface to be a MatrixExpression

template <typename T1, typename T2,
          typename I1 = typename interface<T1>::type,
          typename I2 = typename interface<T2>::type,
          typename Value = typename MatrixAdditionProxy<T1, T2>::value_type>
struct MatrixAdditionInterface;

template <typename T1, typename T2, typename Orient,
          typename I1v, typename I1i,
          typename I2v, typename I2i,
          typename Value>
struct MatrixAdditionInterface<T1, T2,
                             Concepts::DenseMatrix<I1v, Orient, I1i>,
                             Concepts::DenseMatrix<I2v, Orient, I2i>, Value>
{
   typedef Concepts::DenseMatrix<Value, Orient, MatrixAdditionProxy<T1, T2>> type;
   typedef Value value_type;
};

template <typename T1, typename T2,
          typename I1v, typename I1i,
          typename I2v, typename I2i,
          typename Value>
struct MatrixAdditionInterface<T1, T2,
                             Concepts::MatrixExpression<I1v, I1i>,
                             Concepts::MatrixExpression<I2v, I2i>, Value>
{
   typedef Concepts::MatrixExpression<Value, MatrixAdditionProxy<T1, T2>> type;
   typedef Value value_type;
};

template <typename T1, typename T2>
struct interface<MatrixAdditionProxy<T1, T2> >
   : public MatrixAdditionInterface<T1, T2>
{
};

template <typename T1, typename T2>
struct abstract_interface<MatrixAdditionProxy<T1, T2> >
{
   typedef typename matrix_abstract_or<abstract_interface<T1>,
                                       abstract_interface<T2> >::type type;
};

template <typename T1, typename T2>
class MatrixAdditionProxy
{
   public:
   //BOOST_MPL_ASSERT((is_const_proxy_reference_or_immediate<T1>));
   //BOOST_MPL_ASSERT((is_const_proxy_reference_or_immediate<T2>));

      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      typedef T1 reference1;
      typedef T2 reference2;

      typedef typename basic_type<T1>::type matrix1_type;
      typedef typename basic_type<T2>::type matrix2_type;

      // the abstract interface type -
      // shortcut for specializing LinearAlgebra::abstract_interface
   //      typedef typename matrix_abstract_or<abstract_interface<T1>,
   //                                          abstract_interface<T2> >::type abstract_interface;

   //   typedef typename boost::mpl::print<abstract_interface>::type dummy;
   //   typedef typename boost::mpl::print<typename abstract_interface::type>::type dummy2;

      typedef typename interface<T1>::value_type value1_type;
      typedef typename interface<T2>::value_type value2_type;

      typedef Addition<value1_type, value2_type> functor_type;
      typedef typename functor_type::result_type reference;
      typedef typename make_const_reference<reference>::type const_reference;
      typedef typename make_value<reference>::type value_type;

      MatrixAdditionProxy(reference1 x, reference2 y)
         : x_(x), y_(y)
      {
         using LinearAlgebra::size1; using LinearAlgebra::size2;
         DEBUG_PRECONDITION_EQUAL(size1(x), size1(y));
         DEBUG_PRECONDITION_EQUAL(size2(x), size2(y));
      }

      size_type size1() const { return Size1<matrix1_type>()(x_);  }
      size_type size2() const { return Size2<matrix1_type>()(x_); }

   //      const_reference operator[](size_type n) const { return f_(x_[n], y_[n]); }

      reference1 matrix1() const { return x_; }
      reference2 matrix2() const { return y_; }

   private:
      MatrixAdditionProxy& operator=(MatrixAdditionProxy const&); // not implemented

      reference1 x_;
      reference2 y_;
};

// iterators

template <typename T1, typename T2, typename Ti>
struct MatrixAdditionIterate {};

template <typename T1, typename T2>
struct Iterate<MatrixAdditionProxy<T1, T2> >
   : MatrixAdditionIterate<T1, T2, typename interface<MatrixAdditionProxy<T1, T2> >::type> {};

template <typename T1, typename T2, typename Tv, typename TOrient, typename Ti>
struct MatrixAdditionIterate<T1, T2, Concepts::DenseMatrix<Tv, TOrient, Ti>>
{
   typedef typename const_iterator<typename basic_type<T1>::type>::type iter1_type;
   typedef typename const_iterator<typename basic_type<T2>::type>::type iter2_type;
   typedef MatrixAdditionProxy<T1, T2> arg_type;
   typedef typename arg_type::value1_type value1_type;
   typedef typename arg_type::value2_type value2_type;
   typedef MatrixAdditionProxy<T1, T2> const& argument_type;
   typedef MatrixBinaryOuterIterator<iter1_type, iter2_type, typename arg_type::functor_type>
      result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.matrix1()), iterate(x.matrix2()));
   }
};

// binary_transform overloads for addition

template <typename S, typename T, typename F, typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct BinaryTransformMatrix {};

template <typename S, typename T, typename F, typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct BinaryTransformMatrixSemiregular : BinaryTransformMatrix<S, T, F> {};

template <typename S, typename T, typename F, typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct BinaryTransformMatrixRegular : BinaryTransformMatrixSemiregular<S, T, F> {};

// Dispatch to one of the above functors based on F
template <typename S, typename T, typename F,
          bool IsSemiregular = is_semiregular<F>::value,
          bool IsRegular = is_regular<F>::value>
struct BinaryTransformMatrixDispatchFunc : BinaryTransformMatrix<S, T, F> {};

template <typename S, typename T, typename F>
struct BinaryTransformMatrixDispatchFunc<S, T, F, true, true>
   : BinaryTransformMatrixRegular<S, T, F> {};

template <typename S, typename T, typename F>
struct BinaryTransformMatrixDispatchFunc<S, T, F, false, true>
   : BinaryTransformMatrixSemiregular<S, T, F> {};

// a chance to specialize on the functor, before we pass off to the generic version
template <typename S, typename T, typename F, typename Enable = void>
struct BinaryTransformMatrixOverrideFunc : BinaryTransformMatrixDispatchFunc<S, T, F> {};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransform<S, T, F, Concepts::AnyMatrix<Sv, Si>, Concepts::AnyMatrix<Tv, Ti>>
   : BinaryTransformMatrixOverrideFunc<S, T, F> {};

// TODO: we want to commute the arguments here, to get dense matrixs to the left of
// sparse matrixs.
template <typename S, typename T, typename Sv, typename Tv>
struct BinaryTransformMatrixOverrideFunc<S, T, Addition<Sv, Tv> >
{
   typedef typename is_semiregular<Addition<Sv, Tv> >::type semiregular;
   typedef MatrixAdditionProxy<typename make_const_reference<S>::type,
                               typename make_const_reference<T>::type> result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef Addition<Sv, Tv> third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, Addition<Sv, Tv> const& f) const
   { return result_type(x, y); }
};

// a-b is turned into a + (-b)
template <typename S, typename T, typename Sv, typename Tv>
struct BinaryTransformMatrixOverrideFunc<S, T, Subtraction<Sv, Tv> >
{
   typedef typename is_semiregular<Subtraction<Sv, Tv> >::type semiregular;
   typedef typename Negate<T>::result_type NegatedType;
   // this transformation is only valid if NegatedType is a valid
   // proxy-reference, not a temporary.
   //BOOST_MPL_ASSERT((is_proxy_reference_or_immediate<NegatedType>));

   typedef MatrixAdditionProxy<typename make_const_reference<S>::type, NegatedType> result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef Subtraction<Sv, Tv> third_argument_type;

   result_type operator()(S const& x, T const& y) const { return result_type(x, -y); }

   result_type operator()(S const& x, T const& y, Subtraction<Sv, Tv> const& f) const
   { return result_type(x, -y); }
};

// For dense matrices, we ensure that they have the same major index

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransform<S, T, F,
                       Concepts::DenseMatrix<Sv, RowMajor, Si>,
                       Concepts::DenseMatrix<Tv, ColMajor, Ti>>
{
   typedef typename is_semiregular<F>::type semiregular;
   typedef typename is_regular<F>::type regular;
   typedef BinaryTransform<S, typename SwapSortOrder<T>::result_type, F> bt;
   typedef typename bt::result_type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef F third_argument_type;
   result_type operator()(S const& x, T const& y) const
      { return bt()(x, swap_sort_order(y)); }
   result_type operator()(S const& x, T const& y, F f) const
      { return bt()(x, swap_sort_order(y), f); }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransform<S, T, F,
                       Concepts::DenseMatrix<Sv, ColMajor, Si>,
                       Concepts::DenseMatrix<Tv, RowMajor, Ti>>
{
   typedef typename is_semiregular<F>::type semiregular;
   typedef typename is_regular<F>::type regular;
   typedef BinaryTransform<S, typename SwapSortOrder<T>::result_type, F> bt;
   typedef typename bt::result_type result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef F third_argument_type;
   result_type operator()(S const& x, T const& y) const
      { return bt()(x, swap_sort_order(y)); }
   result_type operator()(S const& x, T const& y, F f) const
      { return bt()(x, swap_sort_order(y), f); }
};

// Assign

template <typename LHS, typename T1, typename T2>
struct AssignExpression<LHS&, MatrixAdditionProxy<T1, T2> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixAdditionProxy<T1, T2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x, y.matrix1());
      add(x, y.matrix2());
   }
};

// add

template <typename LHS, typename T1, typename T2>
struct AddExpression<LHS&, MatrixAdditionProxy<T1, T2> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixAdditionProxy<T1, T2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x, y.matrix1());
      add(x, y.matrix2());
   }
};

// subtract

template <typename LHS, typename T1, typename T2>
struct SubtractExpression<LHS&, MatrixAdditionProxy<T1, T2> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef MatrixAdditionProxy<T1, T2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x, y.matrix1());
      subtract(x, y.matrix2());
   }
};

// SwapSortOrder

template <typename T1, typename T2>
struct SwapSortOrder<MatrixAdditionProxy<T1, T2> >
{
   typedef boost::mpl::true_ involutary;
   typedef MatrixAdditionProxy<T1, T2> const& argument_type;
   typedef typename SwapSortOrder<typename basic_type<T1>::type>::result_type r1type;
   typedef typename SwapSortOrder<typename basic_type<T2>::type>::result_type r2type;
   typedef MatrixAdditionProxy<r1type, r2type> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(swap_sort_order(x.matrix1()), swap_sort_order(x.matrix2()));
   }
};


} // namespace LinearAlgebra

#endif
