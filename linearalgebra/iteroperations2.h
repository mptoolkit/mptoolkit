// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/iteroperations2.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

  iteroperations.h

  various algorithms that act on the vector iterator interface,
  some specialized algorithms that act on vectors of vectors (ie. matrices).
  These have the suffix _2, indicating that the operation is applied
  to the nested iterator.

  The default versions look slightly silly; but the point is that
  they can be overloaded and optimized for matrix iterators.

  Algorithms:
  iter_zero_2
  iter_fill_2
  iter_assign_2

  fast_copy
*/

#if !defined(ITEROPERATIONS2_H_HJHUIERH389YU938498PYY98)
#define ITEROPERATIONS2_H_HJHUIERH389YU938498PYY98

#include "iteroperations.h"
#include "matrixptriterator.h"
#include "common/blas1f.h"

namespace LinearAlgebra
{

//
// iter_zero_2
//

template <typename T>
inline
void iter_zero_2(T i)
{
   while (i)
   {
      iter_zero(iterate(*i));
      ++i;
   }
}

//
// iter_fill_2
//

template <typename T, typename V>
inline
void iter_fill_2(T i, V x)
{
   while (i)
   {
      iter_fill(iterate(*i), x);
      ++i;
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
      assign(*dest, *first);
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

// complex<double> type is handled by BLAS
//IMPLEMENT_FAST_COPY_BUILTIN(std::complex<double>)

#undef IMPLEMENT_FAST_COPY_BUILTIN

inline
void fast_copy(double const* first, double const* last, double* dest)
{
   BLAS::dcopy(last-first, first, 1, dest, 1);
}

inline
void fast_copy(std::complex<double> const* first, std::complex<double> const* last,
               std::complex<double>* dest)
{
   BLAS::zcopy(last-first, first, 1, dest, 1);
}

//
// iter_assign_2
//

template <typename I1, typename I2>
inline
void iter_assign_2(I1 i1, I2 i2, vector_iterator_sparse, vector_iterator_sparse)
{
   iter_assign(i1, i2);
}

template <typename I1, typename I2>
inline void
iter_assign_2(I1 const& i1, I2 const& i2)
{
   iter_assign_2(i1, i2, typename I1::category(), typename I2::category());
}

//
// iter_assign_2 specializations for MatrixPtrIterator
//

template <typename Scalar1, typename Orient1, typename Scalar2, typename Orient2>
inline
void iter_assign_2(MatrixPtrIterator<Scalar1, Orient1> const& i1,
                   MatrixPtrIterator<Scalar2, Orient2> const& i2)
{
   DEBUG_PRECONDITION_EQUAL(i1.size1(), i2.size1());
   DEBUG_PRECONDITION_EQUAL(i1.size2(), i2.size2());
   if ((((i1.stride2() == 1 && i1.stride1() == difference_type(i1.size2()))
         && (i2.stride2() == 1 && i2.stride1() == difference_type(i2.size2())))
        || ((i1.stride1() == 1 && i1.stride2() == difference_type(i1.size1()))
            && (i2.stride1() == 1 && i2.stride2() == difference_type(i2.size1()))))
       && (i1.index() == 0) && (i2.index() == 0))
   {
      fast_copy(static_cast<Scalar2 const*>(i2.base()),
                static_cast<Scalar2 const*>(i2.base()) + i1.size1() * i1.size2(),
                i1.base());
   }
   else
   {
      iter_assign(i1, i2);
   }
}

//
// iter_fill_2 specializations for MatrixPtrIterator
//

template <typename T, typename Orient, typename V>
inline
void iter_fill_2(MatrixPtrIterator<T, Orient> i, V x)
{
   // fast path; does the matrix have a linear stride?
   if (i.outer_stride() == i.inner_stride() * difference_type(i.inner_size()))
   {
      iter_fill(VecPtrIterator<T, tagVariable>(i.ptr(), i.size() * i.inner_size(), 0, i.inner_stride()), x);
   }
   // alternatively, is it transposed linear?
   else if (i.inner_stride() == i.outer_stride() * difference_type(i.size()))
   {
      iter_fill(VecPtrIterator<T, tagVariable>(i.ptr(), i.size() * i.inner_size(), 0, i.outer_stride()), x);
   }
   // otherwise do it by hand
   else
   {
      T* I = i.ptr();
      size_type n = 0;
      while (n < i.size())
      {
         iter_fill(VecPtrIterator<T, tagVariable>(I, i.inner_size(), 0, i.inner_stride()), x);
         I += i.outer_stride();
         ++n;
      }
   }
}

//
// iter_max
//

template <typename Iter>
inline
typename Iter::iterator
iter_matrix_max(Iter I)
{
   CHECK(bool(I));
   typedef typename Iter::iterator Inner;
   Inner MaxIter = iter_max(iterate(I));
   ++I;
   while (I)
   {
      Inner This = iter_max(iterate(I));
      if (*MaxIter < *This) MaxIter = This;
      ++I;
   }
   return MaxIter;
}

//
// iter_min
//

template <typename Iter>
inline
typename Iter::iterator
iter_matrix_min(Iter I)
{
   CHECK(bool(I));
   typedef typename Iter::iterator Inner;
   Inner MIter = iter_max(iterate(I));
   ++I;
   while (I)
   {
      Inner This = iter_min(iterate(I));
      if (*This < *MIter) MIter = This;
      ++I;
   }
   return MIter;
}

} // namespace LinearAlgebra

#endif
