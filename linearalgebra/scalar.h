// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/scalar.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
/* -*- C++ -*- $Id$

  scalar.h

  operator specializations for scalar operators
*/

#if !defined(SCALAR_H_DSCHJKLDHVIURGYIURHLUIERHI)
#define SCALAR_H_DSCHJKLDHVIURGYIURHLUIERHI

#include "interface.h"
#include "operations.h"
#include "builtin.h"
#include "common/math_const.h"
#include <complex>
#include <cmath>
#include <cstdlib>

namespace LinearAlgebra
{

//
// AnyScalar
//
// top-level interface for types that act as a scalar
//

template <typename T>
struct AnyScalar
{
   typedef T type;
   typedef T value_type;
};

//
// is_scalar, boolean function to determine if a type
// has interface type of AnyScalar
//

namespace Private
{

template <typename T>
struct is_scalar_helper : boost::mpl::false_ {};

template <typename T>
struct is_scalar_helper<AnyScalar<T> > : boost::mpl::true_ {};

} // namespace Private

template <typename T>
struct is_scalar : Private::is_scalar_helper<typename interface<T>::type> {};

// another misc utility function

template <typename T>
struct is_complex : boost::mpl::false_ { };

template <typename T>
struct is_complex<std::complex<T> > : boost::mpl::true_ { };

//
// For simple cases where the standard arithmetic conventions on commutativity etc hold
// and the value_type is determined through promotion rules, the easiest way to
// enable arithmetic on a given type is to define PromoteTraits for it.
// PromoteTraits is defined for builtin types, but it is explicitly disabled for
// mixed signed arithmetic (a probably misguided attempt to protect the user from
// getting it wrong).  It could be specialized for user types that act like ordinary
// arithmetic types.
//

template <typename T1, typename T2, typename Enable = void>
struct PromoteTraits
{ };

// With gcc 4, numeric_limits<array_type> causes a compile-time failure.  THis causes the
// PromoteTraits specialization to fail.  As a workaround, we have our own version of is_signed,
// which explicitly filters out array types.
template <typename T, typename Enable = void>
struct is_signed : boost::mpl::false_ {};

template <typename T>
struct is_signed<T, typename boost::enable_if<boost::mpl::not_<boost::is_array<T> > >::type>
   : boost::mpl::bool_<std::numeric_limits<T>::is_signed> {};

template <typename T1, typename T2>
struct PromoteTraits<T1, T2, typename boost::enable_if<boost::mpl::and_<
  boost::is_integral<T1>,
  boost::is_integral<T2>,
  boost::mpl::equal_to<is_signed<T1>, is_signed<T2> >
  > >::type>
{
   typedef T1 first_argument_type;
   typedef T2 second_argument_type;
   typedef typename builtin_conversion<T1, T2>::type value_type;
   typedef value_type result_type;
};

template <typename T1, typename T2>
struct PromoteTraits<T1, T2, typename boost::enable_if<boost::mpl::and_<
  boost::is_integral<T1>,
  boost::is_float<T2>
  > >::type>
{
   typedef T1 first_argument_type;
   typedef T2 second_argument_type;
   typedef typename builtin_conversion<T1, T2>::type value_type;
   typedef value_type result_type;
};

template <typename T1, typename T2>
struct PromoteTraits<T1, T2, typename boost::enable_if<boost::mpl::and_<
  boost::is_float<T1>,
  boost::is_arithmetic<T2>
  > >::type>
{
   typedef T1 first_argument_type;
   typedef T2 second_argument_type;
   typedef typename builtin_conversion<T1, T2>::type value_type;
   typedef value_type result_type;
};

template <typename T>
struct PromoteTraits<T, std::complex<T> >
{
   typedef T first_argument_type;
   typedef std::complex<T> second_argument_type;
   typedef std::complex<T> value_type;
   typedef std::complex<T> result_type;
};

template <typename T>
struct PromoteTraits<std::complex<T>, T>
{
   typedef std::complex<T> first_argument_type;
   typedef T second_argument_type;
   typedef std::complex<T> value_type;
   typedef std::complex<T> result_type;
};

template <typename T>
struct PromoteTraits<std::complex<T>, std::complex<T> >
{
   typedef std::complex<T> first_argument_type;
   typedef std::complex<T> second_argument_type;
   typedef std::complex<T> value_type;
   typedef std::complex<T> result_type;
};

//
// Any type that has PromoteTraits defined is assumed to be a scalar
//

template <typename T>
struct interface<T,
                 typename boost::enable_if<
   exists<typename PromoteTraits<T, T>::value_type> >::type>
{
   typedef AnyScalar<T> type;
   typedef T value_type;
};

//
// ScalarProxy
// A proxy to force a type to be a scalar.  For example,
// if x,y are vectors then scalar(x) * y multiplies
// each element of y by x, resulting in a vector of vectors.
//

template <typename T>
struct ScalarProxy
{
   ScalarProxy(T& v) : value(v) {}
   operator T&() { return value; }
   operator T const&() const { return value; }
   T& value;
};

template <typename T>
struct interface<ScalarProxy<T> >
{
   typedef AnyScalar<ScalarProxy<T> > type;
   typedef typename boost::remove_const<T>::type value_type;
};

// Scalar - a functor to coerce an object to a scalar

template <typename T, typename TInterface = typename interface<T>::type>
struct Scalar
{
   typedef T argument_type;
   typedef ScalarProxy<T const> result_type;
   result_type operator()(T const& x) const { return result_type(x); }
};

template <typename T, typename TInterface = typename interface<T>::type>
struct ScalarRef
{
   typedef T& argument_type;
   typedef ScalarProxy<T> result_type;
   result_type operator()(T& x) const { return result_type(x); }
};

template <typename T>
inline
typename Scalar<T>::result_type
scalar(T const& x)
{
   return Scalar<T>()(x);
}

template <typename T>
inline
typename ScalarRef<T>::result_type
scalar(T& x)
{
   return ScalarRef<T>()(x);
}

// Scalar applied to a scalar does nothing

template <typename T>
struct Scalar<T, AnyScalar<T> > : Identity<T> { };

template <typename T>
struct ScalarRef<T, AnyScalar<T> > : IdentityRef<T> { };

template <typename T, typename S>
struct ZeroAllInterface<T&, AnyScalar<S> >
{
   typedef void result_type;
   typedef T& argument_type;
   void operator()(T& v) const
   {
      v = zero<T>();
   }
};

//
// Specializations
//

// by default, transpose on a scalar does nothing.
template <typename T, typename S>
struct TransposeInterface<T, AnyScalar<S> > : Identity<T> {};

// by default, herm on a scalar reduces to the conjugate.
template <typename T, typename S>
struct HermInterface<T, AnyScalar<S> > : Conj<T> {};

template <typename T, typename S>
struct TraceInterface<T, AnyScalar<S> > : Identity<T> {};

template <typename T, typename S>
struct NegateInterface<T, AnyScalar<S>, typename boost::enable_if<
   exists<typename PromoteTraits<T, T>::value_type> >::type>
{
   typedef boost::mpl::true_ stateless;
   typedef boost::mpl::true_ involutary;
   typedef T argument_type;
   typedef T result_type;
   result_type operator()(T const& x) const { return -x; }
};

template <typename T, typename F, typename S>
struct TransformInterface<T, F, AnyScalar<S> >
{
   typedef boost::mpl::true_ stateless;
   typedef T argument_type;
   typedef T first_argument_type;
   typedef F second_argument_type;
   typedef typename std::result_of<F(T)>::type result_type;
   //   typedef typename F::result_type result_type;
   result_type operator()(T const& x, F const& f) const { return f(x); }
   result_type operator()(T const& x) const { return F()(x); }
};

template <typename T, typename F, typename S>
struct TransformInterface<T&, F, AnyScalar<S> >
{
   typedef boost::mpl::true_ stateless;
   typedef T& argument_type;
   typedef T& first_argument_type;
   typedef F second_argument_type;
   typedef typename std::result_of<F(T&)>::type result_type;
   result_type operator()(T& x, F const& f) const { return f(x); }
   result_type operator()(T& x) const { return F()(x); }
};

template <typename S, typename T>
struct AdditionInterface<S, T, AnyScalar<S>, AnyScalar<T>,
                typename boost::enable_if<
   exists<typename PromoteTraits<S, T>::value_type> >::type>
{
   typedef boost::mpl::true_ builtin;
   typedef boost::mpl::true_ semiregular;
   typedef Addition<T, S> commuted_function;

   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x+y; }
};

template <typename S, typename T>
struct SubtractionInterface<S, T, AnyScalar<S>, AnyScalar<T>,
                   typename boost::enable_if<
   exists<typename PromoteTraits<S, T>::value_type> >::type>
{
   typedef boost::mpl::true_ builtin;
   typedef boost::mpl::true_ semiregular;

   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x-y; }
};

template <typename S, typename T>
struct MultiplicationInterface<S, T, AnyScalar<S>, AnyScalar<T>,
                      typename boost::enable_if<
   exists<typename PromoteTraits<S, T>::value_type> >::type>
{
   typedef boost::mpl::true_ builtin;
   typedef Multiplication<T, S> commuted_function;

   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x*y; }
};

template <typename S, typename T, typename Si, typename Ti>
struct DirectProductInterface<S, T, AnyScalar<Si>, AnyScalar<Ti> >
   : MultiplicationInterface<S, T, AnyScalar<Si>, AnyScalar<Ti> > {};

template <typename S, typename T>
struct InnerProdInterface<S, T, AnyScalar<S>, AnyScalar<T>,
                 typename boost::enable_if<
   exists<typename PromoteTraits<S, T>::value_type> >::type>
{
   typedef typename PromoteTraits<S, T>::value_type result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return herm(x)*y; }
};

template <typename S, typename T>
struct EqualToInterface<S, T, AnyScalar<S>, AnyScalar<T>,
               typename boost::enable_if<
   exists<typename PromoteTraits<S, T>::value_type> >::type>
{
   typedef boost::mpl::true_ builtin;

   typedef bool result_type;

   typedef S first_argument_type;
   typedef T second_argument_type;

   result_type operator()(S x, T y) const { return x == y; }
};

// for scalars, assume default-constucted value is zero
template <typename T, typename S>
struct ZeroInterface<T, AnyScalar<S> >
{
   typedef T result_type;
   T operator()() const { return T(); }
};

template <typename T>
struct Zero<T, typename boost::enable_if<
   exists<typename PromoteTraits<T, T>::value_type> >::type>
{
   typedef T result_type;
   T operator()() const { return T(); }
};

template <typename T>
struct IsZero<T, typename boost::enable_if<
   exists<typename PromoteTraits<T, T>::value_type> >::type>
{
   typedef bool result_type;
   bool operator()(T const& x) const { return x == T(); }
};

// by default, the 2-norm and infinity-norm default to the 1-norm.

template <typename T>
struct Norm2<T, AnyScalar<T> >
   : Norm1<T> { };

template <typename T>
struct NormInf<T, AnyScalar<T> >
   : Norm1<T> { };

// the frobenius norm is always the same as the 2-norm here

template <typename T>
struct NormFrobSq<T, AnyScalar<T> > : Norm2Sq<T> {};

//
// specializations for builtins and std::complex
//

// real

template <typename T, typename Res = T, typename Enable = void>
struct ArithmeticIdentity {};

template <typename T, typename Res>
struct ArithmeticIdentity<T, Res, typename boost::enable_if<boost::is_arithmetic<T> >::type>
   : Identity<Res> {};

template <typename T>
struct RealInterface<T, AnyScalar<T> > : ArithmeticIdentity<T> {};

template <typename T>
struct RealInterface<T&, AnyScalar<T> > : ArithmeticIdentity<T, T&> {};

// for std::complex, we do a glorious hack so we can return a reference.
// Hopefully C++0x will make this slightly more legitimate :-)
template <typename T>
struct RealInterface<std::complex<T>, AnyScalar<std::complex<T> > >
{
   typedef T result_type;
   typedef std::complex<T> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return reinterpret_cast<T const*>(&x)[0];
   }
};

template <typename T>
struct RealInterface<std::complex<T>&, AnyScalar<std::complex<T> > >
{
   typedef T& result_type;
   typedef std::complex<T>& argument_type;
   result_type operator()(argument_type x) const
   {
      return reinterpret_cast<T*>(&x)[0];
   }
};

// imag

template <typename T, typename Enable = void>
struct MakeZero {};

template <typename T>
struct MakeZero<T, typename boost::enable_if<is_defined<Zero<T> > >::type>
{
   typedef typename Zero<T>::result_type result_type;
   typedef T argument_type;
   result_type operator()(T const&) const { return Zero<T>()(); }
};

template <typename T>
struct ImagInterface<T, AnyScalar<T>, typename boost::enable_if<boost::is_arithmetic<T> >::type>
   : MakeZero<T> { };

// see comment for Real<std::complex<T> >
template <typename T>
struct Imag<std::complex<T> >
{
   typedef boost::mpl::true_ idempotent;
   typedef T result_type;
   typedef std::complex<T> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return reinterpret_cast<T const*>(&x)[1];
   }
};

template <typename T>
struct Imag<std::complex<T>&>
{
   typedef boost::mpl::true_ idempotent;
   typedef T& result_type;
   typedef std::complex<T>& argument_type;
   result_type operator()(argument_type x) const
   {
      return reinterpret_cast<T*>(&x)[1];
   }
};

// conj

template <typename T>
struct ConjInterface<T, AnyScalar<T>, typename boost::enable_if<boost::is_arithmetic<T> >::type>
   : Identity<T> { };

template <typename T>
struct Conj<std::complex<T> >
{
   typedef boost::mpl::true_ involutary;
   typedef std::complex<T> result_type;
   typedef std::complex<T> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return std::conj(x);
   }
};

// Abs

template <typename T>
struct AbsInterface<T, AnyScalar<T> >
   : Norm2<T> {};

// Arg

template <>
struct ArgInterface<std::complex<float>, AnyScalar<std::complex<float> > >
{
   typedef float result_type;
   typedef std::complex<float> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return std::arg(x);
   }
};

template <>
struct ArgInterface<std::complex<double>, AnyScalar<std::complex<double> > >
{
   typedef double result_type;
   typedef std::complex<double> const& argument_type;
   result_type operator()(argument_type x) const
   {
      return std::arg(x);
   }
};

template <>
struct ArgInterface<float, AnyScalar<float> >
{
   typedef float result_type;
   typedef float argument_type;
   result_type operator()(argument_type x) const
   {
      return x >= 0 ? 0 : float(math_const::pi);
   }
};

template <>
struct ArgInterface<double, AnyScalar<double> >
{
   typedef double result_type;
   typedef double argument_type;
   result_type operator()(argument_type x) const
   {
      return x >= 0 ? 0 : math_const::pi;
   }
};

// Sin

template <>
struct Sin<double>
{
   typedef double result_type;
   typedef double argument_type;
   result_type operator()(argument_type x) const
   {
      return std::sin(x);
   }
};

template <>
struct Sin<std::complex<double> >
{
   typedef std::complex<double> result_type;
   typedef std::complex<double> argument_type;
   result_type operator()(argument_type x) const
   {
      return std::sin(x);
   }
};

// Cos

template <>
struct Cos<double>
{
   typedef double result_type;
   typedef double argument_type;
   result_type operator()(argument_type x) const
   {
      return std::cos(x);
   }
};

template <>
struct Cos<std::complex<double> >
{
   typedef std::complex<double> result_type;
   typedef std::complex<double> argument_type;
   result_type operator()(argument_type x) const
   {
      return std::cos(x);
   }
};


// Exp

template <>
struct Exp<double>
{
   typedef double result_type;
   typedef double argument_type;
   result_type operator()(argument_type x) const
   {
      return std::exp(x);
   }
};

template <>
struct Exp<std::complex<double> >
{
   typedef std::complex<double> result_type;
   typedef std::complex<double> argument_type;
   result_type operator()(argument_type x) const
   {
      return std::exp(x);
   }
};

// Sqrt

template <>
struct Sqrt<double>
{
   typedef double result_type;
   typedef double argument_type;
   result_type operator()(argument_type x) const
   {
      return std::sqrt(x);
   }
};

// No need for transpose here, by default it does nothing for scalar types.

// No need for herm here, defaults to conj for scalar types.

//
// 1-norm and squared 2-norm
//

template <typename T>
struct Norm2Sq<T, AnyScalar<T>, typename boost::enable_if<boost::is_arithmetic<T> >::type>
{
   typedef double result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
   {
      return x * x;
   }
};

template <typename T>
struct Norm2Sq<std::complex<T> >
{
   typedef T result_type;
   typedef std::complex<T> argument_type;
   result_type operator()(argument_type const& x) const
   {
      return std::norm(x);
   }
};

template <typename T>
struct Norm1<T, AnyScalar<T>, typename boost::enable_if<boost::is_arithmetic<T> >::type>
{
   typedef T result_type;
   typedef T argument_type;
   result_type operator()(T const& x) const
   {
      return std::abs(x);
   }
};

template <typename T>
struct Norm1<std::complex<T> >
{
   typedef T result_type;
   typedef std::complex<T> argument_type;
   result_type operator()(argument_type const& x) const
   {
      // FIXME: this isn't as accurate as it could be; use pythag() instead
      return std::sqrt(std::norm(x));
   }
};

// equal

template <typename T, typename U, typename TolType>
struct EqualInterface<T, U, TolType, AnyScalar<T>, AnyScalar<U> >
{
   typedef bool result_type;
   typedef T first_argument_type;
   typedef U second_argument_type;

   EqualInterface(TolType Tol = default_tolerance()) : Tol_(Tol) {}

   bool operator()(T const& x, U const& y) const
   {
      return norm_inf(x-y)
         <= Tol_ + 2 * std::numeric_limits<TolType>::epsilon() * (norm_inf(x) + norm_inf(y));
   }

   private:
      TolType Tol_;
};

// assign for scalars defaults to operator=

template <typename S, typename T>
struct AssignInterface<S&, T, AnyScalar<S>, AnyScalar<T> >
{
   typedef void result_type;
   typedef S& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      x = y;
   }
};

template <typename S, typename T>
struct AddInterface<S&, T, AnyScalar<S>, AnyScalar<T> >
{
   typedef void result_type;
   typedef S& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      x += y;
   }
};

template <typename S, typename T>
struct SubtractInterface<S&, T, AnyScalar<S>, AnyScalar<T> >
{
   typedef void result_type;
   typedef S& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      x -= y;
   }
};

template <typename S, typename T>
struct MultiplyInterface<S&, T, AnyScalar<S>, AnyScalar<T> >
{
   typedef void result_type;
   typedef S& first_argument_type;
   typedef T const& second_argument_type;
   result_type operator()(S& x, T const& y) const
   {
      x *= y;
   }
};

} // namespace LinearAlgebra

#endif
