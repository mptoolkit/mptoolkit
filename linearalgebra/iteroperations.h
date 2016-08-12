// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/iteroperations.h
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
  iteroperations.h

  various algorithms that act on the vector iterator interface.

  Algorithms:
  iter_zero
  iter_fill
  iter_assign
  iter_add
  iter_subtract
  iter_equal_to
  iter_norm_1
  iter_norm_2_sq
  iter_norm_frob_sq
  iter_norm_inf
  iter_inner_prod
  iter_inner_prod_dense
  iter_coefficient_inner_prod
  iter_coefficient_inner_prod_dense
  iter_min
  iter_max
  iter_sum

  iter_inner_prod_sparse is defined in another header
*/

#if !defined(MPTOOLKIT_LINEARALGEBRA_ITEROPERATIONS_H)
#define MPTOOLKIT_LINEARALGEBRA_ITEROPERATIONS_H

#include "scalar.h"
#include "vectorinterface.h"
#include "vecptriterator.h"
#include <cstring>

namespace LinearAlgebra
{

//
// iter_zero
//

template <typename T>
inline
void
iter_zero(T i)
{
   while (i)
   {
      zero_all(*i);
      ++i;
   }
}

// Assuming we have IEEE arithmetic, the representation of 0.0 is all zeros.
inline
void iter_zero(VecPtrIterator<float> i)
{
   std::memset(i.ptr(), 0, i.size() * sizeof(float));
}

inline
void iter_zero(VecPtrIterator<double> i)
{
   std::memset(i.ptr(), 0, i.size() * sizeof(double));
}

//
// iter_fill
//

template <typename I, typename V>
inline
void iter_fill(I i, V x)
{
   while (i)
   {
      *i = x;
      ++i;
   }
}

// for pointer types, we can use std::fill which may be faster
// (or may not - surprisingly STL implementations don't seem to 
//  do these optimizations)
template <typename T, typename V>
inline
void iter_fill(VecPtrIterator<T> i, V x)
{
   std::fill(i.ptr(), i.end(), x);
}

// even better, use memset for char types
template <typename V>
inline
void iter_fill(VecPtrIterator<char> i, V x)
{
   std::memset(i.ptr(), (unsigned int)(char(x)), i.size());
}

template <typename V>
inline
void iter_fill(VecPtrIterator<unsigned char> i, V x)
{
   std::memset(i.ptr(), (unsigned char)(x), i.size());
}

template <typename V>
inline
void iter_fill(VecPtrIterator<signed char> i, V x)
{
   std::memset(i.ptr(), (unsigned int)((signed char)(x)), i.size());
}

template <typename V>
void iter_fill(VecPtrIterator<float> i, V x)
{
   float f = x;
   if (f == 0)
   {
      std::memset(i.ptr(), 0, i.size() * sizeof(float));
   }
   else
   {
      std::fill(i.ptr(), i.end(), f);
   }
}

template <typename V>
void iter_fill(VecPtrIterator<double> i, V x)
{
   double f = x;
   if (f == 0)
   {
      std::memset(i.ptr(), 0, i.size() * sizeof(double));
   }
   else
   {
      std::fill(i.ptr(), i.end(), f);
   }
}

template <typename V>
void iter_fill(VecPtrIterator<int> i, V x)
{
   int f = x;
   if (f == 0)
   {
      std::memset(i.ptr(), 0, i.size() * sizeof(int));
   }
   else
   {
      std::fill(i.ptr(), i.end(), f);
   }
}

template <typename V>
void iter_fill(VecPtrIterator<unsigned int> i, V x)
{
   unsigned int f = x;
   if (f == 0)
   {
      std::memset(i.ptr(), 0, i.size() * sizeof(unsigned int));
   }
   else
   {
      std::fill(i.ptr(), i.end(), f);
   }
}

// TODO: iter_fill for a stride iterator

//
// iter_assign
//

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_sparse)
{
   iter_zero(i1);
   while (i2)
   {
      add(i1[i2.index()], *i2);
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_injective)
{
   iter_zero(i1);
   while (i2)
   {
      assign(i1[i2.index()], *i2);
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_ordered)
{
   // hmm, this override may be slower than the injective version in practice
   while (i2)
   {
      size_type i2Index = i2.index();
      while (i1.index() < i2Index)
      {
         assign(*i1, zero<typename I1::value_type>());
	 ++i1;
      }
      assign(*i1, *i2);
      ++i2;
      ++i1;
   }
   while (i1)
   {
      assign(*i1, zero<typename I1::value_type>());
      ++i1;
   }
}

template <typename I1, typename I2>
void iter_assign(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_dense)
{
   while (i2)
   {
      assign(*i1, *i2);
      ++i1;
      ++i2;
   }
}

template <typename I1, typename I2>
inline
void iter_assign(I1 const& i1, I2 const& i2)
{
   iter_assign(i1, i2, typename I1::category(), typename I2::category());
}

//
// iter_add
//

template <typename I1, typename I2>
void iter_add(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_sparse)
{
   while (i2)
   {
      add(i1[i1.index()], *i2);
      ++i1;
   }
}

template <typename I1, typename I2>
void iter_add(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_dense)
{
   while (i1)
   {
      DEBUG_CHECK(bool(i2));
      add(*i1, *i2);
      ++i1;
      ++i2;
   }
   DEBUG_CHECK(!i2);
}

template <typename I1, typename I2>
void iter_add(I1 const& i1, I2 const& i2)
{
   iter_add(i1, i2, typename I1::category(), typename I2::category());
}

//
// iter_subtract
//

template <typename I1, typename I2>
void iter_subtract(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_sparse)
{
   while (i2)
   {
      subtract(i1[i1.index()], *i2);
      ++i1;
   }
}

template <typename I1, typename I2>
void iter_subtract(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_dense)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   while (i1)
   {
      subtract(*i1, *i2);
      ++i1;
      ++i2;
   }
}

template <typename I1, typename I2>
void iter_subtract(I1 const& i1, I2 const& i2)
{
   iter_subtract(i1, i2, typename I1::category(), typename I2::category());
}

//
// iter_equal_to
//

template <typename I1, typename I2>
bool iter_equal_to(I1 i1, I2 i2, vector_iterator_dense, vector_iterator_dense)
{
   if (i1 && i2)
   {
      if (i1.index() != i2.index()) return false;
      while (i1 && i2)
      {
	 if (*i1 != *i2) return false;
	 ++i1;
	 ++i2;
      }
   }
   return !(i1 || i2);
}

template <typename I1, typename I2>
bool iter_equal_to(I1 i1, I2 i2, vector_iterator_ordered, vector_iterator_ordered)
{
   while (i1 && i2)
   {
      if (i1.index() != i2.index()) return false;
      if (*i1 != *i2) return false;
      ++i1;
      ++i2;
   }
   return !(i1 || i2);
}

template <typename I1, typename I2>
bool iter_equal_to(I1 i1, I2 i2, vector_iterator_injective, vector_iterator_hashed)
{
   if (i1.size() != i2.size()) return false;

   size_type i1Size = 0;
   while (i1)
   {
      if (*i1 != i2(i1.index())) return false;
      ++i1;
      ++i1Size;
   }
   return i1Size == i2.size();
   return true;
}

template <typename I1, typename I2>
bool iter_equal_to(I1 i1, I2 i2, vector_iterator_hashed, vector_iterator_injective)
{
   size_type i2Size = 0;
   while (i2)
   {
      if (i1(i2.index()) != *i2) return false;
      ++i2;
      ++i2Size;
   }
   return i2Size == i1.size();
   return true;
}

template <typename I1, typename I2>
bool iter_equal_to(I1 i1, I2 i2, vector_iterator_hashed, vector_iterator_hashed)
{
   // TODO: there is probably a faster way to implement this

   I1 I = i1;
   while (I)
   {
      if (*I != i2(I.index())) return false;
      ++I;
   }

   while (i2)
   {
      if (*i2 != i1(i2.index())) return false;
      ++i2;
   }

   return true;
}

template <typename I1, typename I2>
inline
bool iter_equal_to(I1 const& i1, I2 const& i2)
{
   return iter_equal_to(i1, i2, typename I1::category(), typename I2::category());
}

//
// iter_norm_1
//

template <typename I>
typename Norm1<typename I::value_type>::result_type 
iter_norm_1(I i, vector_iterator_injective)
{
   typedef typename Norm1<typename I::value_type>::result_type result_type;
   if (!i) return zero_or_die<result_type>();

   result_type x = norm_1(*i);
   ++i;
   while (i)
   {
      x += norm_1(*i);
      ++i;
   }
   return x;
}

template <typename I>
inline
typename Norm1<typename I::value_type>::result_type 
iter_norm_1(I const& i)
{
   return iter_norm_1(i, typename I::category());
}

//
// iter_norm_2_sq
//

template <typename I>
typename Norm2Sq<typename I::value_type>::result_type 
iter_norm_2_sq(I i, vector_iterator_injective)
{
   typedef typename Norm2Sq<typename I::value_type>::result_type result_type;
   if (!i) return zero_or_die<result_type>();

   result_type x = norm_2_sq(*i);
   ++i;
   while (i)
   {
      x += norm_2_sq(*i);
      ++i;
   }
   return x;
}

template <typename I>
inline
typename Norm2Sq<typename I::value_type>::result_type 
iter_norm_2_sq(I const& i)
{
   return iter_norm_2_sq(i, typename I::category());
}

//
// iter_norm_2
//

template <typename I>
inline
typename Norm2Sq<typename I::value_type>::result_type 
iter_norm_2(I const& i)
{
   using std::sqrt;
   return sqrt(iter_norm_2_sq(i, typename I::category()));
}


//
// iter_norm_frob_sq
//

template <typename I>
typename NormFrobSq<typename I::value_type>::result_type 
iter_norm_frob_sq(I i, vector_iterator_injective)
{
   typedef typename NormFrobSq<typename I::value_type>::result_type result_type;
   if (!i) return zero_or_die<result_type>();

   result_type x = norm_frob_sq(*i);
   ++i;
   while (i)
   {
      x += norm_frob_sq(*i);
      ++i;
   }
   return x;
}

template <typename I>
inline
typename NormFrobSq<typename I::value_type>::result_type 
iter_norm_frob_sq(I const& i)
{
   return iter_norm_frob_sq(i, typename I::category());
}

//
// iter_norm_frob
//

template <typename I>
inline
typename NormFrobSq<typename I::value_type>::result_type 
iter_norm_frob(I const& i)
{
   using std::sqrt;
   return sqrt(iter_norm_frob_sq(i, typename I::category()));
}


//
// iter_norm_inf
//
// we do a bit of a hack to avoid some square roots when evaluating
// the norm_inf of a complex vector.  Since we know that norm_inf(complex<T> x)
// is the same as norm_2(x), we can use norm_2_sq(x) instead and take
// the square root only at the end.

template <typename I>
typename boost::disable_if<is_complex<typename I::value_type>,
			   typename NormInf<typename I::value_type>::result_type>::type
iter_norm_inf(I i, vector_iterator_injective)
{
   using std::max;
   typedef typename NormInf<typename I::value_type>::result_type result_type;
   if (!i) return zero_or_die<result_type>();

   result_type Max = norm_inf(*i);
   ++i;
   while (i)
   {
      Max = max(Max, norm_inf(*i));
      ++i;
   }
   return Max;
}

template <typename I>
typename boost::enable_if<is_complex<typename I::value_type>,
			  typename NormInf<typename I::value_type>::result_type>::type
iter_norm_inf(I i, vector_iterator_injective)
{
   using std::max;
   typedef typename NormInf<typename I::value_type>::result_type result_type;
   if (!i) return zero_or_die<result_type>();

   result_type Max = norm_2_sq(*i);
   ++i;
   while (i)
   {
      Max = max(Max, norm_2_sq(*i));
      ++i;
   }
   return std::sqrt(Max);
}

template <typename I>
inline
typename NormInf<typename I::value_type>::result_type 
iter_norm_inf(I const& i)
{
   return iter_norm_inf(i, typename I::category());
}

//
// iter_inner_prod
//

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_dense, vector_iterator_dense)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      DEBUG_PRECONDITION(!bool(i2));
      return zero_or_die<value_type>();
   }

   DEBUG_PRECONDITION(bool(i2));
   value_type Result(f(*i1, *i2));
   ++i1; ++i2;
   while (i1)
   {
      DEBUG_PRECONDITION(bool(i2));
      add(Result, f(*i1, *i2));
      ++i1; ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_dense, vector_iterator_sparse)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1[i2.index()], *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1[i2.index()], *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_sparse, vector_iterator_dense)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(*i1, i2[i1.index()]));
   ++i1;
   while (i1)
   {
      add(Result, f(*i1, i2[i1.index()]));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_dense, vector_iterator_hashed)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1[i2.index()], *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1[i2.index()], *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_hashed, vector_iterator_dense)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(*i1, i2[i1.index()]));
   ++i1;
   while (i1)
   {
      add(Result, f(*i1, i2[i1.index()]));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_dense, vector_iterator_ordered)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1[i2.index()], *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1[i2.index()], *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_ordered, vector_iterator_dense)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(*i1, i2[i1.index()]));
   ++i1;
   while (i1)
   {
      add(Result, f(*i1, i2[i1.index()]));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_ordered, vector_iterator_sparse)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_sparse, vector_iterator_ordered)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(*i1, i2(i1.index())));
   ++i1;
   while (i1)
   {
      add(Result, f(*i1, i2(i1.index())));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_hashed, vector_iterator_sparse)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_sparse, vector_iterator_hashed)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(*i1, i2(i1.index())));
   ++i1;
   while (i1)
   {
      add(Result, f(*i1, i2(i1.index())));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_hashed, vector_iterator_ordered)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_ordered, vector_iterator_hashed)
{
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(*i1, i2(i1.index())));
   ++i1;
   while (i1)
   {
      add(Result, f(*i1, i2(i1.index())));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_hashed, vector_iterator_hashed)
{
   // TODO: this could be made more efficient for some cases
   // by choosing which iterator to increment vs lookup; ie
   // it would be better to increment the iterator that has the fewest elements.

   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      add(Result, f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 i1, I2 i2, Func f, vector_iterator_ordered, vector_iterator_ordered)
{
   typedef typename make_value<typename Func::result_type>::type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1 || !i2) return zero_or_die<value_type>();

   while (i1.index() != i2.index())
   {
      if (i1.index() < i2.index())
      {
         ++i1;
         if (!i1) return zero_or_die<value_type>();
      }
      else
      {
         ++i2;
         if (!i2) return zero_or_die<value_type>();
      }
   }

   DEBUG_CHECK_EQUAL(i1.index(), i2.index());
   result_type Result(f(*i1, *i2));
   ++i1; ++i2;

   while (i1 && i2)
   {
      if (i1.index() < i2.index())
      {
         ++i1;
         while (i1 && i1.index() < i2.index())
            ++i1;
      }
      else
      {
         ++i2;
         while (i2 && i2.index() < i1.index())
            ++i2;
      }

      if (i1 && i2)
      {
         DEBUG_CHECK_EQUAL(i1.index(), i2.index());
         Result += f(*i1, *i2);
         ++i1; ++i2;
      }
   }
   return Result;
}

template <typename I1, typename I2, typename Func>
inline
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod(I1 const& i1, I2 const& i2, Func const& f)
{
   return iter_inner_prod(i1, i2, f, typename I1::category(), typename I2::category());
}

//
// iter_coefficient_inner_prod
//

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func const& f,
                            vector_iterator_dense, vector_iterator_dense)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      DEBUG_PRECONDITION(!bool(i2));
      return zero_or_die<value_type>();
   }

   DEBUG_PRECONDITION(bool(i2));
   value_type Result(cf(i1.index()) * f(*i1, *i2));
   ++i1; ++i2;
   while (i1)
   {
      DEBUG_PRECONDITION(bool(i2));
      {
         result_type x = cf(i1.index());
         if (x != 0)
            //      if (result_type x = cf(i1.index()) != 0)
         add(Result, x * f(*i1, *i2));
      }
      //      add(Result, cf(i1.index()) * f(*i1, *i2));
      ++i1; ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_dense, vector_iterator_sparse)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1[i2.index()], *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
      if (!is_zero(x))
         add(Result, x * f(i1[i2.index()], *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf,
                            Func f, vector_iterator_sparse, vector_iterator_dense)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i1.index()) * f(*i1, i2[i1.index()]));
   ++i1;
   while (i1)
   {
      cf_result x(cf(i1.index()));
         add(Result, x * f(*i1, i2[i1.index()]));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_dense, vector_iterator_hashed)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1[i2.index()], *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
         add(Result, x * f(i1[i2.index()], *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf,
                            Func f, vector_iterator_hashed, vector_iterator_dense)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i1.index()) * f(*i1, i2[i1.index()]));
   ++i1;
   while (i1)
   {
      cf_result x(cf(i1.index()));
         add(Result, x * f(*i1, i2[i1.index()]));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_dense, vector_iterator_ordered)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1[i2.index()], *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
         add(Result, x * f(i1[i2.index()], *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf,
                            Func f, vector_iterator_ordered, vector_iterator_dense)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i1.index()) * f(*i1, i2[i1.index()]));
   ++i1;
   while (i1)
   {
      cf_result x(cf(i1.index()));
      if (!is_zero(x))
         add(Result, x * f(*i1, i2[i1.index()]));

      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_ordered, vector_iterator_sparse)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
      if (!is_zero(x))
         add(Result, x * f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf,
                            Func f, vector_iterator_sparse, vector_iterator_ordered)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i1.index()) * f(*i1, i2(i1.index())));
   ++i1;
   while (i1)
   {
      cf_result x(cf(i1.index()));
      if (!is_zero(x))
         add(Result, x * f(*i1, i2(i1.index())));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_hashed, vector_iterator_sparse)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
      if (!is_zero(x))
         add(Result, x * f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_sparse, vector_iterator_hashed)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result( cf(i1.index()) * f(*i1, i2(i1.index())));
   ++i1;
   while (i1)
   {
      cf_result x(cf(i1.index()));
      if (!is_zero(x))
         add(Result,  x * f(*i1, i2(i1.index())));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_hashed, vector_iterator_ordered)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
      if (!is_zero(x))
         add(Result, x * f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, Func f, 
                            vector_iterator_ordered, vector_iterator_hashed)
{
   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1)
   {
      return zero_or_die<value_type>();
   }

   value_type Result( cf(i1.index()) * f(*i1, i2(i1.index())));
   ++i1;
   while (i1)
   {
      cf_result x(cf(i1.index()));
      if (!is_zero(x))
         add(Result,  x * f(*i1, i2(i1.index())));
      ++i1;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, 
                            Func f, vector_iterator_hashed, vector_iterator_hashed)
{
   // TODO: this could be made more efficient for some cases
   // by choosing which iterator to increment vs lookup; ie
   // it would be better to increment the iterator that has the fewest elements.

   typedef typename CF::result_type cf_result;
   typedef typename Func::result_type result_type;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i2)
   {
      return zero_or_die<value_type>();
   }

   value_type Result(cf(i2.index()) * f(i1(i2.index()), *i2));
   ++i2;
   while (i2)
   {
      cf_result x(cf(i2.index()));
      if (!is_zero(x))
         add(Result, x * f(i1(i2.index()), *i2));
      ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 i1, I2 i2, CF const& cf, 
                            Func const& f, vector_iterator_ordered, vector_iterator_ordered)
{
   typedef typename result_value<Func>::type result_type;
   typedef typename CF::result_type cf_result;
   typedef typename make_value_with_zero<result_type>::type value_type;
   if (!i1 || !i2) return zero_or_die<value_type>();

   while (i1.index() != i2.index())
   {
      if (i1.index() < i2.index())
      {
         ++i1;
         if (!i1) return zero_or_die<value_type>();
      }
      else
      {
         ++i2;
         if (!i2) return zero_or_die<value_type>();
      }
   }

   DEBUG_CHECK_EQUAL(i1.index(), i2.index());

   cf_result x(cf(i2.index()));
   while (is_zero(x))
   {
      ++i1; ++i2;
      if (!i1 || !i2) return value_type();

      while (i1.index() != i2.index())
      {
         if (i1.index() < i2.index())
         {
            ++i1;
            if (!i1) return value_type();
         }
         else
         {
            ++i2;
            if (!i2) return value_type();
         }
      }
      x = cf(i2.index());
   }

   value_type Result(x * f(*i1, *i2));
   ++i1; ++i2;

   while (i1 && i2)
   {
      while (i1.index() != i2.index())
      {
         if (i1.index() < i2.index())
         {
            ++i1;
            if (!i1) return Result;
         }
         else
         {
            ++i2;
            if (!i2) return Result;
         }
      }

      DEBUG_CHECK_EQUAL(i1.index(), i2.index());
      cf_result x(cf(i2.index()));
      if (!is_zero(x))
         Result += x * f(*i1, *i2);

      ++i1; ++i2;
   }
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
inline
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod(I1 const& i1, I2 const& i2, CF const& cf, Func const& f)
{
   return iter_coefficient_inner_prod(i1, i2, cf, f, 
                                      typename I1::category(), typename I2::category());
}

//
// iter_max
//

template <typename Iter>
inline
Iter iter_max(Iter I, vector_iterator_injective)
{
   CHECK(bool(I));
   Iter MaxIter = I;
   ++I;
   while (I)
   {
      if (*MaxIter < *I) MaxIter = I;
      ++I;
   }
   return MaxIter;
}

template <typename Iter>
inline
Iter iter_max(Iter const& I)
{
   return iter_max(I, typename Iter::category());
}

//
// iter_min
//

template <typename Iter>
inline
Iter iter_min(Iter I, vector_iterator_injective)
{
   CHECK(bool(I));
   Iter MinIter = I;
   ++I;
   while (I)
   {
      if (*I < *MinIter) MinIter = I;
      ++I;
   }
   return MinIter;
}

template <typename Iter>
inline
Iter iter_min(Iter const& I)
{
   return iter_min(I, typename Iter::category());
}

//
// iter_sum
//

template <typename I>
typename make_value_with_zero<typename I::value_type>::type 
iter_sum(I i, vector_iterator_injective)
{
   typedef typename make_value_with_zero<typename I::value_type>::type result_type;
   if (!i) return zero_or_die<result_type>();

   result_type x = *i;
   ++i;
   while (i)
   {
      x += *i;
      ++i;
   }
   return x;
}

template <typename I>
inline
typename make_value_with_zero<typename I::value_type>::type 
iter_sum(I const& i)
{
   return iter_sum(i, typename I::category());
}

} // namespace LinearAlgebra

//#define NO_BLAS


#if !defined(NO_BLAS)
#include "iteroperations_blas.h"
#endif

#endif
