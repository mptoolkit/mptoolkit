// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/functors.h
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

// Generic function objects for wrapping operators and functions.
// These are useful for nested operations, eg a nested tensor product can
// be evaluated as tensor_prod(A, B, TensorProd())
// versus the default tensor_prod(A, B, Multiplication())

#if !defined(MPTOOLKIT_BLAS_FUNCTORS_H)
#define MPTOOLKIT_BLAS_FUNCTORS_H

#include "blas/number_traits.h"

namespace blas
{

struct Multiplication
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return x*y;
   }
};

struct Addition
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return x+y;
   }
};

struct Subtraction
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return x-y;
   }
};

struct InnerProd
{
   template <typename T, typename U>
   auto operator()(T const& x, U const& y) const
   {
      return inner_prod_nested(x,y, *this);
   }
   // 3-argument version takes the result as a parameter
   template <typename T, typename U, typename V>
   auto operator()(T const& x, U const& y, V&& Result) const
   {
      return inner_prod_nested(x,y,Result, *this);
   }
   template <typename T, typename U, typename V>
   auto add(T const& x, U const& y, V&& Result) const
   {
      add_inner_prod_nested(x,y, Result, *this);
   }
};
 
struct InnerProdNested
{
   template <typename T, typename U, typename Nested>
   auto operator()(T const& x, U const& y, Nested&& Nest) const
   {
      return inner_prod_nested(x,y, std::forward<Nested>(Nest));
   }
   // 3-argument version takes the result as a parameter
   template <typename T, typename U, typename V, typename Nested>
   auto operator()(T const& x, U const& y, V&& Result, Nested&& Nest) const
   {
      return inner_prod_nested(x,y,Result,std::forward<Nested>(Nest));
   }
   template <typename T, typename U, typename V, typename Nested>
   auto add(T const& x, U const& y, V&& Result, Nested&& Nest) const
   {
      add_inner_prod_nested(x,y, Result, std::forward<Nested>(Nest));
   }
};

template <typename T, typename U>
inline
auto inner_prod(T const& x, U const& y)
{
   return inner_prod_nested(x,y,InnerProd());
}

template <typename T, typename U, typename R>
inline
auto inner_prod(T const& x, U const& y, R&& Result)
{
   return inner_prod_nested(x,y,std::forward<R>(Result), InnerProd());
}

template <typename T, typename U, typename R>
inline
auto add_inner_prod(T const& x, U const& y, R&& Result)
{
   return add_inner_prod_nested(x,y,std::forward<R>(Result), InnerProd());
}

//
// inner_prod for floating point and complex types
//

template <typename T, typename Nested>
inline
std::enable_if_t<blas::is_numeric_v<T>, T>
inner_prod_nested(T x, T y, Nested&&)
{
   using std::conj;
   return conj(x)*y;
}

template <typename T, typename Nested>
inline
std::enable_if_t<blas::is_numeric_v<T>, void>
inner_prod_nested(T x, T y, T& r, Nested&&)
{
   using std::conj;
   r = conj(x)*y;
}

template <typename T, typename Nested>
inline
std::enable_if_t<blas::is_numeric_v<T>, void>
add_inner_prod_nested(T x, T y, T& r, Nested&&)
{
   using std::conj;
   r += conj(x)*y;
}

// Frobenius norm

template <typename T>
inline
std::enable_if_t<std::is_floating_point<T>::value, T>
norm_frob(T x)
{
   return std::abs(x);
}

template <typename T>
inline
T
norm_frob(std::complex<T> const& x)
{
   return std::hypot(x.real(), x.imag());
}

template <typename T>
inline
std::enable_if_t<std::is_floating_point<T>::value, T>
norm_frob_sq(T x)
{
   return x*x;
}

template <typename T>
inline
T
norm_frob_sq(std::complex<T> const& x)
{
   return std::norm(x);
}

// trace

template <typename T>
inline
std::enable_if_t<blas::is_numeric<T>::value, T>
trace(T x)
{
   return x;
}

// conj

template <typename T>
inline
std::enable_if_t<blas::is_real_v<T>, void>
inplace_conj(T& x)
{
}

template <typename T>
inline
std::enable_if_t<blas::is_complex_v<T>, void>
inplace_conj(T& x)
{
   using std::conj;
   x = conj(x);
}

// herm

template <typename T>
inline
std::enable_if_t<blas::is_real_v<T>, T>
herm(T const& x)
{
   using std::conj;
   return x;
}

template <typename T>
inline
std::enable_if_t<blas::is_complex_v<T>, T>
herm(T const& x)
{
   using std::conj;
   return conj(x);
}

} // namespace blas

using blas::norm_frob;
using blas::norm_frob_sq;
using blas::inner_prod;
using blas::inplace_conj;
using blas::herm;

#endif
