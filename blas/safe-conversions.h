// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/safe-conversions.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
// 'Safe' conversions of scalar types.
// Declares the metafunction convert<U>(T const& x) that
// converts a value of type T to a value of type U.
// This is a subset of possible implicit conversions, plus some
// other conversions that are not automatic in C++ (for example,
// conversions from integral types to complex<double>).  We disallow
// conversions that may be dangerous in some circumstances, for example
// conversion float <-> double can be dangerous in either direction, depending
// on the circumstances.  (Converting float to double can be dangerous if it
// originates in some constant value that should be obtained with higher precision,
// eg double(sqrtf(2.0)) is erroneous.)
//

#if !defined(MPTOOLKIT_BLAS_SAFE_CONVERSIONS_H)
#define MPTOOLKIT_BLAS_SAFE_CONVERSIONS_H

#include "config.h"
#include <type_traits>

namespace blas
{

namespace detail
{

// dummy template that is always false so that
// the static assert in convert_impl depends on the template parameters
template <typename T, typename U>
struct fail : public std::false_type {};

// conversion from T to U
template <typename T, typename U>
struct convert_impl
{
   static_assert(fail<T,U>::value, "conversion is not allowed");
};

template <typename T>
struct convert_impl<T,T>
{
   static T convert(T const& x) { return x; }
};

#define DECLARE_CONVERSION(a,b) template <> struct convert_impl<a,b> \
   { static b convert(a const& x) { return (b)(x); } }

#define DECLARE_CONVERSION_INTERMEDIATE(a,b,i) template <> struct convert_impl<a,b> \
   { static b convert(a const& x) { return (b)(x); } }

DECLARE_CONVERSION(int,long);
DECLARE_CONVERSION(int,long long);
DECLARE_CONVERSION(long,long long);
DECLARE_CONVERSION(int,float);
DECLARE_CONVERSION(int,double);
DECLARE_CONVERSION(long,double);
DECLARE_CONVERSION(double,std::complex<double>);
DECLARE_CONVERSION_INTERMEDIATE(int,std::complex<double>,double);
DECLARE_CONVERSION_INTERMEDIATE(long,std::complex<double>,double);

#if defined(HAVE_FLOAT128)
DECLARE_CONVERSION(int,float128);
DECLARE_CONVERSION(long,float128);
DECLARE_CONVERSION(long long,float128);
DECLARE_CONVERSION_INTERMEDIATE(int,std::complex<float128>, float128);
DECLARE_CONVERSION_INTERMEDIATE(long,std::complex<float128>, float128);
DECLARE_CONVERSION_INTERMEDIATE(long long,std::complex<float128>, float128);
DECLARE_CONVERSION(float128,std::complex<float128>);
#endif

#undef DECLARE_CONVERSION
#undef DECLARE_CONVERSION_INTERMEDIATE

} // namespace detail

template <typename T, typename U, typename Enable = void>
struct is_safe_conversion : public std::false_type {};

template <typename T, typename U>
struct is_safe_conversion<T, U, std::enable_if<
   std::is_same< U,
                 decltype(detail::convert_impl<T,U>::convert(std::declval<T>()))>::value>>
   : public std::true_type {};

template <typename T, typename U>
using is_safe_conversion_v = typename is_safe_conversion<T,U>::value;

template <typename T, typename U>
inline
T safe_convert(U const& x)
{
   return detail::convert_impl<U,T>::convert(x);
}

} // namespace blas

#endif
