// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectoraddition.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  vectoraddition.h

  A proxy class & iterators for binary expressions of vectors

  Created 2005-01-10 Ian McCulloch
*/

#if !defined(VECTORADDITION_H_KDSHJCKYH48Y98)
#define VECTORADDITION_H_KDSHJCKYH48Y98

#include "vectoroperationsbase.h"
#include "vectorbinaryiterator.h"
#include <boost/static_assert.hpp>

namespace LinearAlgebra
{

template <typename T1, typename T2>
class VectorAdditionProxy
{
   public:
      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      // the abstract interface type -
      // shortcut for specializing LinearAlgebra::abstract_interface
      typedef vector_abstract_or<typename abstract_interface<T1>::type,
                                 typename abstract_interface<T2>::type> abstract_interface;

      typedef typename interface<T1>::value_type value1_type;
      typedef typename interface<T2>::value_type value2_type;

      typedef Addition<value1_type, value2_type> functor_type;
      typedef typename functor_type::result_type reference;
      typedef typename make_const_reference<reference>::type const_reference;
      typedef typename make_value<reference>::type value_type;
      typedef typename make_const_reference<T1>::type reference1;
      typedef typename make_const_reference<T2>::type reference2;

      VectorAdditionProxy(reference1 x, reference2 y)
         : x_(x), y_(y) {} //{ CHECK_EQUAL(Size<T1>()(x_), Size<T2>()(y_)); }

      size_type size() const { using LinearAlgebra::size; return size(x_); }
   // return Size<T1>()(x_); }

      const_reference operator[](size_type n) const { return f_(x_[n], y_[n]); }

      reference1 vector1() const { return x_; }
      reference2 vector2() const { return y_; }

   private:
      reference1 x_;
      reference2 y_;
};

// interface

template <typename T1, typename T2,
          typename I1 = typename interface<T1>::type,
          typename I2 = typename interface<T2>::type,
          typename Value = typename VectorAdditionProxy<T1, T2>::value_type>
struct VectorAdditionInterface;

template <typename T1, typename T2,
          typename I1v, typename I1i,
          typename I2v, typename I2i,
          typename Value>
struct VectorAdditionInterface<T1, T2,
                             DENSE_VECTOR(I1v, I1i),
                             DENSE_VECTOR(I2v, I2i), Value>
{
   typedef DENSE_VECTOR(Value, void) type;
   typedef Value value_type;
};

template <typename T1, typename T2,
          typename I1v, typename I1i,
          typename I2v, typename I2i,
          typename Value>
struct VectorAdditionInterface<T1, T2,
                             VECTOR_EXPRESSION(I1v, I1i),
                             VECTOR_EXPRESSION(I2v, I2i), Value>
{
   typedef VECTOR_EXPRESSION(Value, void) type;
   typedef Value value_type;
};

template <typename T1, typename T2>
struct interface<VectorAdditionProxy<T1, T2> >
   : public VectorAdditionInterface<T1, T2>
{
};

// iterators

template <typename T1, typename T2>
struct Iterate<VectorAdditionProxy<T1, T2> >
{
   typedef typename const_iterator<T1>::type iter1_type;
   typedef typename const_iterator<T2>::type iter2_type;
   typedef VectorAdditionProxy<T1, T2> const& argument_type;
   typedef typename VectorAdditionProxy<T1, T2>::functor_type func_type;
   typedef VectorBinaryIterator<iter1_type, iter2_type, func_type> result_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.vector1()), iterate(x.vector2()));
   }
};

// binary_transform overloads for addition

template <typename S, typename T, typename F, typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct BinaryTransformVector {};

template <typename S, typename T, typename F,  typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct BinaryTransformVectorSemiregular : BinaryTransformVector<S, T, F> {};

template <typename S, typename T, typename F,  typename Si = typename interface<S>::type,
          typename Ti = typename interface<T>::type>
struct BinaryTransformVectorRegular : BinaryTransformVectorSemiregular<S, T, F> {};

// Dispatch to one of the above functors based on F
template <typename S, typename T, typename F,
          bool IsSemiregular = is_semiregular<F>::value,
          bool IsRegular = is_regular<F>::value>
struct BinaryTransformVectorDispatchFunc : BinaryTransformVector<S, T, F> {};

template <typename S, typename T, typename F>
struct BinaryTransformVectorDispatchFunc<S, T, F, true, true>
   : BinaryTransformVectorRegular<S, T, F> {};

template <typename S, typename T, typename F>
struct BinaryTransformVectorDispatchFunc<S, T, F, false, true>
   : BinaryTransformVectorSemiregular<S, T, F> {};

// a chance to specialize on the functor, before we pass off to the generic version
template <typename S, typename T, typename F, typename Enable = void>
struct BinaryTransformVectorOverrideFunc : BinaryTransformVectorDispatchFunc<S, T, F> {};

template <typename S, typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct BinaryTransform<S, T, F, ANY_VECTOR(Sv, Si), ANY_VECTOR(Tv, Ti) >
   : BinaryTransformVectorOverrideFunc<S, T, F> {};

// TODO: we want to commute the arguments here, to get dense vectors to the left of
// sparse vectors.
template <typename S, typename T, typename Sv, typename Tv>
struct BinaryTransformVectorOverrideFunc<S, T, Addition<Sv, Tv> >
{
   typedef boost::mpl::true_ semiregular;
   typedef VectorAdditionProxy<S, T> result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef Addition<Sv, Tv> third_argument_type;
   result_type operator()(S const& x, T const& y) const { return result_type(x, y); }

   result_type operator()(S const& x, T const& y, Addition<Sv, Tv> const& f) const
   { return result_type(x, y); }
};

// a-b is turned into a + (-b)
template <typename S, typename T, typename Sv, typename Tv>
struct BinaryTransformVectorOverrideFunc<S, T, Subtraction<Sv, Tv> >
{
   typedef boost::mpl::true_ semiregular;
   typedef typename Negate<T>::result_type NegatedType;
   typedef VectorAdditionProxy<S, NegatedType> result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef Subtraction<Sv, Tv> third_argument_type;

   // this transformation is only valid if NegatedType is a valid
   // proxy-reference, not a temporary.
   BOOST_STATIC_ASSERT(is_proxy_reference<NegatedType>::value);

   result_type operator()(S const& x, T const& y) const { return result_type(x, -y); }

   result_type operator()(S const& x, T const& y, Subtraction<Sv, Tv> const& f) const
   { return result_type(x, -y); }
};

// Assign

template <typename LHS, typename T1, typename T2>
struct AssignExpression<LHS&, VectorAdditionProxy<T1, T2> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorAdditionProxy<T1, T2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      assign(x, y.vector1());
      add(x, y.vector2());
   }
};

// add

template <typename LHS, typename T1, typename T2>
struct AddExpression<LHS&, VectorAdditionProxy<T1, T2> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorAdditionProxy<T1, T2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      add(x, y.vector1());
      add(x, y.vector2());
   }
};

// subtract

template <typename LHS, typename T1, typename T2>
struct SubtractExpression<LHS, VectorAdditionProxy<T1, T2> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef VectorAdditionProxy<T1, T2> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      subtract(x, y.vector1());
      subtract(x, y.vector2());
   }
};

} // namespace LinearAlgebra

#endif
