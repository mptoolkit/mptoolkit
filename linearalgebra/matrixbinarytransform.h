// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixbinarytransform.h
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
  matrixbinarytransform.h

  A proxy class & iterators for binary expressions of matrices

  Created 2005-03-15 Ian McCulloch
*/

#if !defined(MATRIXBINARYTRANSFORM_H_KDSHJCKYH48Y98)
#define MATRIXBINARYTRANSFORM_H_KDSHJCKYH48Y98

#include "matrixaddition.h"

namespace LinearAlgebra
{

template <typename T1, typename T2, typename F>
class MatrixBinaryTransformProxy;

template <typename T1, typename T2, typename F>
struct abstract_interface<MatrixBinaryTransformProxy<T1, T2, F> >
{
   typedef typename matrix_abstract_or<abstract_interface<T1>,
                                       abstract_interface<T2> >::type type;
};

template <typename T1, typename T2, typename F>
class MatrixBinaryTransformProxy
{
   public:
   //BOOST_MPL_ASSERT((is_const_proxy_reference_or_immediate<T1>));
   //BOOST_MPL_ASSERT((is_const_proxy_reference_or_immediate<T2>));

      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      typedef T1 reference1;
      typedef T2 reference2;
      typedef F functor_type;

      typedef typename basic_type<T1>::type matrix1_type;
      typedef typename basic_type<T2>::type matrix2_type;

      typedef typename interface<T1>::value_type value1_type;
      typedef typename interface<T2>::value_type value2_type;

      typedef typename functor_type::result_type reference;
      typedef typename make_const_reference<reference>::type const_reference;
      typedef typename make_value<reference>::type value_type;

      MatrixBinaryTransformProxy(reference1 x, reference2 y, functor_type f = functor_type())
         : x_(x), y_(y), f_(f)
      {
         using LinearAlgebra::size1; using LinearAlgebra::size2;
         DEBUG_PRECONDITION_EQUAL(size1(x), size1(y));
         DEBUG_PRECONDITION_EQUAL(size2(x), size2(y));
      }

      size_type size1() const { return Size1<matrix1_type>()(x_);  }
      size_type size2() const { return Size2<matrix1_type>()(x_); }

      reference1 matrix1() const { return x_; }
      reference2 matrix2() const { return y_; }

      functor_type const& functor() const { return f_; }

   private:
      MatrixBinaryTransformProxy& operator=(MatrixBinaryTransformProxy const&); // not implemented

      reference1 x_;
      reference2 y_;
      functor_type f_;
};

// interface

template <typename T1, typename T2, typename F,
          typename I1 = typename interface<T1>::type,
          typename I2 = typename interface<T2>::type,
          typename Value = typename MatrixBinaryTransformProxy<T1, T2, F>::value_type>
struct MatrixBinaryTransformInterface;

template <typename T1, typename T2, typename F, typename Orient,
          typename I1v, typename I1i,
          typename I2v, typename I2i,
          typename Value>
struct MatrixBinaryTransformInterface<T1, T2, F,
                                      Concepts::DenseMatrix<I1v, Orient, I1i>,
                                      Concepts::DenseMatrix<I2v, Orient, I2i>, Value>
{
   typedef Concepts::DenseMatrix<Value, Orient, void> type;
   typedef Value value_type;
};

template <typename T1, typename T2, typename F>
struct interface<MatrixBinaryTransformProxy<T1, T2, F> >
   : public MatrixBinaryTransformInterface<T1, T2, F>
{
};

// iterators

template <typename T1, typename T2, typename F, typename Ti>
struct MatrixBinaryTransformIterate {};

template <typename T1, typename T2, typename F>
struct Iterate<MatrixBinaryTransformProxy<T1, T2, F> >
   : MatrixBinaryTransformIterate<T1, T2, F,
                                  typename interface<MatrixBinaryTransformProxy<T1, T2, F> >::type> {};

template <typename T1, typename T2, typename F, typename Tv, typename TOrient, typename Ti>
struct MatrixBinaryTransformIterate<T1, T2, F, Concepts::DenseMatrix<Tv, TOrient, Ti>>
{
   typedef typename const_iterator<typename basic_type<T1>::type>::type iter1_type;
   typedef typename const_iterator<typename basic_type<T2>::type>::type iter2_type;
   typedef MatrixBinaryTransformProxy<T1, T2, F> arg_type;
   typedef typename arg_type::value1_type value1_type;
   typedef typename arg_type::value2_type value2_type;
   typedef MatrixBinaryTransformProxy<T1, T2, F> const& argument_type;
   typedef MatrixBinaryOuterIterator<iter1_type, iter2_type, F>
      result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.matrix1()), iterate(x.matrix2()), x.functor());
   }
};

// BinaryTransform specializations

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti, typename Orient>
struct BinaryTransformMatrix<S, T, F,
                             Concepts::DenseMatrix<Sv, Orient, Si>,
                             Concepts::DenseMatrix<Tv, Orient, Ti>>
{
   typedef MatrixBinaryTransformProxy<typename make_const_reference<S>::type,
                               typename make_const_reference<T>::type, F> result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef F third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, F const& f) const
   { return result_type(x, y, f); }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti, typename Orient>
struct BinaryTransformMatrixSemiregular<S, T, F,
                                        Concepts::CompressedOuterMatrix<Sv, Orient, Si>,
                                        Concepts::CompressedOuterMatrix<Tv, Orient, Ti>>
{
   typedef MatrixBinaryTransformProxy<typename make_const_reference<S>::type,
                               typename make_const_reference<T>::type, F> result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef F third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, F const& f) const
   { return result_type(x, y, f); }
};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransformMatrixSemiregular<S, T, F,
                                        Concepts::DiagonalMatrix<Sv,Si>,
                                        Concepts::DiagonalMatrix<Tv,Ti>>
{
   typedef MatrixBinaryTransformProxy<typename make_const_reference<S>::type,
                               typename make_const_reference<T>::type, F> result_type;
   typedef S const& first_argument_type;
   typedef T const& second_argument_type;
   typedef F third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, F const& f) const
   { return result_type(x, y, f); }
};

// SwapSortOrder

template <typename T1, typename T2, typename F>
struct SwapSortOrder<MatrixBinaryTransformProxy<T1, T2, F> >
{
   typedef boost::mpl::true_ involutary;
   typedef MatrixBinaryTransformProxy<T1, T2, F> const& argument_type;
   typedef typename SwapSortOrder<typename basic_type<T1>::type>::result_type r1type;
   typedef typename SwapSortOrder<typename basic_type<T2>::type>::result_type r2type;
   typedef MatrixBinaryTransformProxy<r1type, r2type, F> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(swap_sort_order(x.matrix1()), swap_sort_order(x.matrix2()), x.functor());
   }
};


} // namespace LinearAlgebra

#endif
