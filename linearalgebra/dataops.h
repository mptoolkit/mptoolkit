// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/dataops.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
  dataops.h

  various algorithms that act on a reduced iterator interface.  The emphasis is on speed.

  The restricted iterator interface is:
  value_type defined
  pre-increment (returns reference to *this)
  operator* (convertible to value_type)
  operator==
  operator!=

  This also includes some basic iterator classes, for which optimized versions
  of the algorithms are possible.

  This could make good use of typeof()

  Algorithms:
  fast_fill
  fast_copy
  fast_copy_scaled
  fast_copy_transformed
  fast_add
  fast_add_scaled
  fast_add_transformed
  fast_subtract
  fast_subtract_scaled
  fast_subtract_transformed
  fast_multiply
  fast_equal
  fast_norm_2_sq
  fast_norm_2
  fast_norm_inf
  fast_min
  fast_max
  fast_inner_prod
  fast_scalar_prod

  dot
  norm2sq
  norm2
  conj (supplements the std:: version)
*/

#if !defined(DATAOPS_H_CHJFRH43RERJ334TY89ERVHUI83REY89)
#define DATAOPS_H_CHJFRH43RERJ334TY89ERVHUI83REY89

#include <iterator>
#include <complex>
#include <algorithm>
#include "operations.h"
#include "iteratortypes.h"
#include "common/trace.h"
//#include "matrixiterators.h"
#include "common/numerics.h"

namespace ops
{

using namespace LinearAlgebra;
using namespace numerics;

//
// fast_fill
//

template <typename Iter, typename Scalar>
void fast_fill(Iter first, Iter last, Scalar s)
{
   while (first != last)
   {
      *first = s;
      ++first;
   }
}

// for pointer types, we can use std::fill which may be faster
// (or may not - surprisingly most STL implementations don't seem to
//  do these optimizations)
template <typename T, typename Scalar>
inline
void fast_fill(T* first, T* last, Scalar s)
{
   std::fill(first, last, s);
}

// even better, use memset for char types
inline
void fast_fill(char* first, char* last, char c)
{
   memset(first, c, last-first);
}

inline
void fast_fill(unsigned char* first, unsigned char* last, unsigned char c)
{
   memset(first, c, last-first);
}

inline
void fast_fill(signed char* first, signed char* last, signed char c)
{
   memset(first, c, last-first);
}

// Assuming we have IEEE arithmetic, the representation of 0.0 is all zeros.
// This is an important enough special case to optimize for it.
inline
void fast_fill(float* first, float* last, float s)
{
   if (s == 0)
   {
      memset(first, 0, (last-first) * sizeof(float));
   }
   else std::fill(first, last, s);
}

inline
void fast_fill(double* first, double* last, double s)
{
   if (s == 0)
   {
      memset(first, 0, (last-first) * sizeof(double));
   }
   else std::fill(first, last, s);
}

// fast_fill for a StrideIterator
template <typename Iter, typename Scalar>
void fast_fill(StrideIterator<Iter> first, StrideIterator<Iter> last, Scalar s)
{
   DEBUG_CHECK_EQUAL(first.stride(), last.stride());
   if (first.stride() == 1)
   {
      fast_fill(first.base(), last.base(), s);
   }
   else
   {
      while (first != last)
      {
         *first = s;
         ++first;
      }
   }
}

//
// fast_copy
//

template <typename Iter1, typename Iter2>
void fast_copy(Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest = *first;
      ++dest; ++first;
   }
}

// for pointer types we can use std::copy
template <typename T1, typename T2>
inline
void fast_copy(T1* first, T1* last, T2* dest)
{
   std::copy(first, last, dest);
}

// Most implementations of std::copy fail to specialize for builtins,
// so we do it ourselves.
// if we had pod_traits we could generalize this
#define IMPLEMENT_FAST_COPY_BUILTIN(BuiltinT)                                   \
inline                                                                          \
void fast_copy(BuiltinT const* first, BuiltinT const* last, BuiltinT* dest)     \
{                                                                               \
   memcpy(dest, first, (last-first) * sizeof(BuiltinT));                        \
}

IMPLEMENT_FAST_COPY_BUILTIN(char)
IMPLEMENT_FAST_COPY_BUILTIN(signed char)
IMPLEMENT_FAST_COPY_BUILTIN(unsigned char)
IMPLEMENT_FAST_COPY_BUILTIN(int)
IMPLEMENT_FAST_COPY_BUILTIN(unsigned int)
IMPLEMENT_FAST_COPY_BUILTIN(long)
IMPLEMENT_FAST_COPY_BUILTIN(unsigned long)
#if defined(USE_LONG_LONG)
IMPLEMENT_FAST_COPY_BUILTIN(long long)
IMPLEMENT_FAST_COPY_BUILTIN(unsigned long long)
#endif

// double type is handled by BLAS
//IMPLEMENT_FAST_COPY_BUILTIN(double)

// storage format for std::complex<double> is not yet mandated to be equivalent
// to double[2], but will be in next standard
IMPLEMENT_FAST_COPY_BUILTIN(std::complex<double>)

#undef IMPLEMENT_FAST_COPY_BUILTIN

// fast_copy on a TransformIterator turns into fast_copy_transform
template <typename BaseIter, typename Functor, typename Iter2>
inline
void fast_copy(TransformIterator<BaseIter, Functor> const& first,
               TransformIterator<BaseIter, Functor> const& last,
               Iter2 const& dest)
{
   fast_copy_transform(first.func(), first.base(), last.base(), dest);
}

//
// fast_copy_scaled
//

// Pass x by value to avoid aliasing problems.
template <typename Scalar, typename Iter1, typename Iter2>
void fast_copy_scaled(Scalar x, Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest = x * (*first);
      ++dest; ++first;
   }
}

//
// fast_copy_transform
//

// Pass x by value to avoid aliasing problems.
template <typename Functor, typename Iter1, typename Iter2>
void fast_copy_transform(Functor const& f, Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest = f(*first);
      ++dest; ++first;
   }
}

#if 0
// fast_copy_transform on a functor that is
// BinderFirst<BinaryOperator<Multiplication, T1, T2> >
// gets transformed into fast_copy_scaled
template <typename T1, typename T2, typename Iter1, typename Iter2>
inline
void fast_copy_transform(BinderFirst<BinaryOperator<Multiplication, T1, T2> > const& f,
                         Iter1 const& first, Iter1 const& last,
                         Iter2 const& dest)
{
   fast_copy_scaled(f.b, first, last, dest);
}
#endif

//
// fast_add
//

template <typename Iter1, typename Iter2>
void fast_add(Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest += *first;
      ++dest; ++first;
   }
}

// fast_add on a TransformIterator turns into fast_add_transform
template <typename BaseIter, typename Functor, typename Iter2>
inline
void fast_add(TransformIterator<BaseIter, Functor> const& first,
              TransformIterator<BaseIter, Functor> const& last,
              Iter2 const& dest)
{
   fast_add_transform(first.func(), first.base(), last.base(), dest);
}

//
// fast_add_scaled
//

// Pass x by value to avoid aliasing problems.
template <typename Scalar, typename Iter1, typename Iter2>
void fast_add_scaled(Scalar x, Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest += x * (*first);
      ++dest; ++first;
   }
}

//
// fast_add_transform
//

template <typename Functor, typename Iter1, typename Iter2>
void fast_add_transform(Functor f, Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest += f(*first);
      ++dest; ++first;
   }
}

#if 0
// fast_add_transform on a functor that is
// BinderFirst<BinaryOperator<Multiplication, T1, T2> >
// gets transformed into fast_add_scaled
template <typename T1, typename T2, typename Iter1, typename Iter2>
inline
void fast_add_transform(BinderFirst<BinaryOperator<Multiplication, T1, T2> > const& f,
                        Iter1 const& first, Iter1 const& last,
                        Iter2 const& dest)
{
   fast_add_scaled(f.b, first, last, dest);
}
#endif

// fast_add_transform on a functor that is Negate<F> gets
// transformed into fast_subtract.
// There is an obscure corner case with these transformations: if type T is
// different from Iter1::value_type as well as Iter2::value_type, then
// this alters the sequence of conversions.  However, it would be unusual
// (and probably a bug) if T was different from Iter1::value_type.
// (eg. consider Iter1::value_type == Iter2::value_type == double, and T = int !)
template <typename T, typename Iter1, typename Iter2>
inline
void fast_add_transform(Negate<T>, Iter1 const& first, Iter1 const& last,
                        Iter2 const& dest)
{
   fast_subtract(first, last, dest);
}

template <typename T, typename F, typename Iter1, typename Iter2>
inline
void fast_add_transform(UnaryComposer<Negate<T>, F> const& f, Iter1 const& first,
                        Iter1 const& last,
                        Iter2 const& dest)
{
   fast_subtract_transform(f.second, first, last, dest);
}

//
// fast_substract
//

template <typename Iter1, typename Iter2>
void fast_subtract(Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest -= *first;
      ++dest; ++first;
   }
}

//
// fast_subtract_scaled
//

template <typename Scalar, typename Iter1, typename Iter2>
void fast_subtract_scaled(Scalar x, Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest -= x * (*first);
      ++dest; ++first;
   }
}

//
// fast_subtract_transform
//

template <typename Functor, typename Iter1, typename Iter2>
void fast_subtract_transform(Functor f, Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      *dest -= f(*first);
      ++dest; ++first;
   }
}

#if 0
// fast_subtract_transform on a functor that is
// BinderFirst<BinaryOperator<Multiplication, T1, T2> >
// gets transformed into fast_subtract_scaled
template <typename T1, typename T2, typename Iter1, typename Iter2>
inline
void fast_subtract_transform(BinderFirst<BinaryOperator<Multiplication, T1, T2> > const& f,
                             Iter1 const& first, Iter1 const& last,
                             Iter2 const& dest)
{
   fast_subtract_scaled(f.b, first, last, dest);
}
#endif

// fast_subtract_transform on a functor that is Negate<F> gets
// transformed into fast_add
template <typename T, typename Iter1, typename Iter2>
inline
void fast_subtract_transform(Negate<T>, Iter1 const& first, Iter1 const& last,
                             Iter2 const& dest)
{
   fast_add(first, last, dest);
}

template <typename T, typename F, typename Iter1, typename Iter2>
inline
void fast_subtract_transform(UnaryComposer<Negate<T>, F> const& f, Iter1 const& first,
                             Iter1 const& last,
                             Iter2 const& dest)
{
   fast_add_transform(f.second, first, last, dest);
}

//
// fast_multiply
//

// multiply by scalar
template <typename Iter1, typename Scalar>
void fast_multiply(Iter1 first, Iter1 last, Scalar x)
{
   while (first != last)
   {
      *first *= x;
      ++first;
   }
}

// equal using pre-increment only.
template <typename Iter1, typename Iter2>
bool fast_equal(Iter1 first, Iter1 last, Iter2 dest)
{
   while (first != last)
   {
      if (*dest != *first) return false;
      ++dest; ++first;
   }
   return true;
}

//
// conj
//

#define IMPLEMENT_BUILTIN_CONJ(T)               \
inline                                          \
T conj(T x)                                     \
{                                               \
  return x;                                     \
}

IMPLEMENT_BUILTIN_CONJ(char)
IMPLEMENT_BUILTIN_CONJ(signed char)
IMPLEMENT_BUILTIN_CONJ(unsigned char)
IMPLEMENT_BUILTIN_CONJ(short)
IMPLEMENT_BUILTIN_CONJ(unsigned short)
IMPLEMENT_BUILTIN_CONJ(int)
IMPLEMENT_BUILTIN_CONJ(unsigned int)
IMPLEMENT_BUILTIN_CONJ(long)
IMPLEMENT_BUILTIN_CONJ(unsigned long)
#if defined(USE_LONG_LONG)
IMPLEMENT_BUILTIN_CONJ(long long)
IMPLEMENT_BUILTIN_CONJ(unsigned long long)
#endif
IMPLEMENT_BUILTIN_CONJ(double)

#undef IMPLEMENT_BUILTIN_CONJ

inline
std::complex<double> conj(std::complex<double> x)
{
   return std::conj(x);
}

//
// fast_2norm2 is the square of the 2-norm.  It always returns a 'double'
//

template <typename Iter>
double fast_norm_2_sq(Iter first, Iter last)
{
   double acc = 0;
   while (first != last)
   {
      acc += norm_2_sq(*first);
      ++first;
   }
   return acc;
}

//
// fast_norm2.  This is defined as the square root of fast_norm2sq,
// it should usually be unnecessary to overload/specialize.
//

template <typename Iter>
inline
double fast_norm_2(Iter first, Iter last)
{
   return std::sqrt(fast_norm_2_sq(first, last));
}

//
// fast_norm_inf
//

template <typename Iter>
inline
double fast_norm_inf(Iter first, Iter last)
{
   double acc = 0;
   while (first != last)
   {
      acc = std::max<double>(acc, norm_inf(*first));
      ++first;
   }
   return acc;
}

//
// fast_min
//

template <typename AccType, typename Iter>
inline
AccType fast_min(AccType acc, Iter first, Iter last)
{
   while (first != last)
   {
      acc = std::min(acc, *first);
      ++first;
   }
   return acc;
}

//
// fast_max
//

template <typename AccType, typename Iter>
inline
AccType fast_max(AccType acc, Iter first, Iter last)
{
   while (first != last)
   {
      acc = std::max(acc, *first);
      ++first;
   }
   return acc;
}

//
// fast_inner_prod
//

template <typename ResultType, typename Iter1, typename Iter2>
void fast_inner_prod(ResultType& x, Iter1 first1, Iter1 last1, Iter2 first2)
{
   while (first1 != last1)
   {
      x += (*first1) * (*first2);
      ++first1;
      ++first2;
   }
}

//
// fast_scalar_prod
//

template <typename ResultType, typename Iter1, typename Iter2>
void fast_scalar_prod(ResultType& x, Iter1 first1, Iter1 last1, Iter2 first2)
{
   while (first1 != last1)
   {
      x += scalar_prod(*first1, *first2);
      ++first1;
      ++first2;
   }
}

} // namespace ops

// bring in the BLAS specializations
#if !defined(NO_BLAS)
#include "dataops_blas.h"
#endif

#endif
