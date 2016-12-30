// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/dataops_blas.h
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

  dataops_blas.h

  Optimizations for the algorithms in dataops.h using BLAS level 1.

  There is a small combinatorial explosion here since we need to define versions for
  all combinations of StrideIterator and normal iterators.  Unfortunately there isn't
  much we can do about this; although one possibility would be to use recursive templates
  in the iterator classes themselves - new class PointerIterator would wrap a T* and inherit from
  StrideIterator, overriding the stride to be 1.
*/

#if !defined(DATAOPS_BLAS_H_DSC843Y8HFE8H8OHJOHFDSIUFLRHWLO)
#define DATAOPS_BLAS_H_DSC843Y8HFE8H8OHJOHFDSIUFLRHWLO

#include "common/blas1f.h"

namespace ops
{

using LinearAlgebra::StrideIterator;

// For BLAS, the required pointer is to the lowest element in memory.  For negative stride,
// this is actually the 'back' of the array.  array_base() is a helper function to
// determine this.
template <typename T>
inline
T* array_base(StrideIterator<T*> first, difference_type n)
{
   return first.stride() < 0 ? first.base() + (n-1)*first.stride() : first.base();
}

// this version instead takes a [first, last) pair
template <typename T>
inline
T* array_base(StrideIterator<T*> first, StrideIterator<T*> last)
{
   DEBUG_PRECONDITION(first.stride() == last.stride());
   return first.stride() < 0 ? last.base()-first.stride() : first.base();
}

// no fast_fill equivalent in BLAS

//
// fast_copy
//

inline
void fast_copy(double const* first, double const* last, double* dest)
{
   BLAS::dcopy(last-first, first, 1, dest, 1);
}

inline
void fast_copy(double const* first, double const* last, StrideIterator<double*> dest)
{
   BLAS::dcopy(last-first, first, 1, array_base(dest, last-first), dest.stride());
}

inline
void fast_copy(StrideIterator<double const*> first,
               StrideIterator<double const*> last,
               double* dest)
{
   BLAS::dcopy(last-first, array_base(first, last), first.stride(), dest, 1);
}

inline
void fast_copy(StrideIterator<double const*> first,
               StrideIterator<double const*> last,
               StrideIterator<double*> dest)
{
   difference_type n = last-first;
   BLAS::dcopy(n, array_base(first, last), first.stride(),
               array_base(dest, n), dest.stride());
}

// there is no scaled dcopy in BLAS - we could implement fast_copy_scaled with
// dcopy followed by dscal, but its probably not worth it.

//
// fast_add
//

inline
void fast_add(double const* first, double const* last, double* dest)
{
   BLAS::daxpy(last-first, 1.0, first, 1, dest, 1);
}

inline
void fast_add(StrideIterator<double const*> first,
              StrideIterator<double const*> last, double* dest)
{
   BLAS::daxpy(last-first, 1.0, array_base(first, last), first.stride(), dest, 1);
}

inline
void fast_add(double const* first, double const* last, StrideIterator<double*> dest)
{
   BLAS::daxpy(last-first, 1.0, first, 1, array_base(dest, last-first), dest.stride());
}

inline
void fast_add(StrideIterator<double const*> first, StrideIterator<double const*> last,
              StrideIterator<double*> dest)
{
   difference_type n = last-first;
   BLAS::daxpy(n, 1.0, array_base(first, last), first.stride(),
               array_base(dest, n), dest.stride());
}

//
// fast_add_scaled
//

inline
void fast_add_scaled(double x, double const* first, double const* last, double* dest)
{
   BLAS::daxpy(last-first, x, first, 1, dest, 1);
}

inline
void fast_add_scaled(double x,
                     StrideIterator<double const*> first,
                     StrideIterator<double const*> last,
                     double* dest)
{
   BLAS::daxpy(last-first, x, array_base(first, last), first.stride(), dest, 1);
}

inline
void fast_add_scaled(double x, double const* first, double const* last,
                     StrideIterator<double*> dest)
{
   BLAS::daxpy(last-first, x, first, 1, array_base(dest, last-first), dest.stride());
}

inline
void fast_add_scaled(double x,
                     StrideIterator<double const*> first,
                     StrideIterator<double const*> last,
                     StrideIterator<double*> dest)
{
   difference_type n = last-first;
   BLAS::daxpy(n, x, array_base(first, last), first.stride(),
               array_base(dest, n), dest.stride());
}

//
// fast_subtract
//

inline
void fast_subtract(double const* first, double const* last, double* dest)
{
   BLAS::daxpy(last-first, -1.0, first, 1, dest, 1);
}

inline
void fast_subtract_scaled(double x, double const* first, double const* last, double* dest)
{
   BLAS::daxpy(last-first, -x, first, 1, dest, 1);
}

//
// fast_multiply
//

inline
void fast_multiply(double* first, double* last, double x)
{
   BLAS::dscal(last-first, x, first, 1);
}

inline
void fast_multiply(StrideIterator<double*> first, StrideIterator<double*> last, double x)
{
   BLAS::dscal(last-first, x, array_base(first, last), first.stride());
}

inline
double fast_dot(double const* first, double const* last, double const* other)
{
   return BLAS::ddot(last-first, first, 1, other, 1);
}

//
// fast_norm_2_sq
//

inline
double fast_norm_2_sq(double const* first, double const* last)
{
   double x = BLAS::dnrm2(last-first, first, 1);
   return x*x;
}

inline
double fast_norm_2_sq(StrideIterator<double const*> first,
                      StrideIterator<double const*> last)
{
   double x = BLAS::dnrm2(last-first, array_base(first, last), first.stride());
   return x*x;
}

//
// fast_norm_2
//

inline
double fast_norm_2(double const* first, double const* last)
{
   return BLAS::dnrm2(last-first, first, 1);
}

inline
double fast_norm_2(StrideIterator<double const*> first, StrideIterator<double const*> last)
{
   return BLAS::dnrm2(last-first, array_base(first, last), first.stride());
}

//
// fast_inner_prod
//

inline
void fast_inner_prod(double& x, double const* first1, double const* last1, double const* first2)
{
   x += BLAS::ddot(last1-first1, first1, 1, first2, 1);
}

inline
void fast_inner_prod(double& x, double const* first1, double const* last1,
                     StrideIterator<double const*> first2)
{
   difference_type n = last1-first1;
   x += BLAS::ddot(n, first1, 1, array_base(first2, n), first2.stride());
}

inline
void fast_inner_prod(double& x, StrideIterator<double const*> first1,
                     StrideIterator<double const*> last1,
                     double const* first2)
{
   x += BLAS::ddot(last1-first1, array_base(first1, last1), first1.stride(), first2, 1);
}

inline
void fast_inner_prod(double& x, StrideIterator<double const*> first1,
                     StrideIterator<double const*> last1,
                     StrideIterator<double const*> first2)
{
   difference_type n = last1-first1;
   x += BLAS::ddot(n, array_base(first1, n), first1.stride(), array_base(first2, n), first2.stride());
}

//
// fast_scalar_prod
//

inline
void fast_scalar_prod(double& x, double const* first1, double const* last1, double const* first2)
{
   x += BLAS::ddot(last1-first1, first1, 1, first2, 1);
}

inline
void fast_scalar_prod(double& x, double const* first1, double const* last1,
                     StrideIterator<double const*> first2)
{
   difference_type n = last1-first1;
   x += BLAS::ddot(n, first1, 1, array_base(first2, n), first2.stride());
}

inline
void fast_scalar_prod(double& x, StrideIterator<double const*> first1,
                     StrideIterator<double const*> last1,
                     double const* first2)
{
   x += BLAS::ddot(last1-first1, array_base(first1, last1), first1.stride(), first2, 1);
}

inline
void fast_scalar_prod(double& x, StrideIterator<double const*> first1,
                     StrideIterator<double const*> last1,
                     StrideIterator<double const*> first2)
{
   difference_type n = last1-first1;
   x += BLAS::ddot(n, array_base(first1, n), first1.stride(), array_base(first2, n), first2.stride());
}

} // namespace ops

#endif
