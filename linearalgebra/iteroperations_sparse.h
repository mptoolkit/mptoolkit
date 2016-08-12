// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/iteroperations_sparse.h
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

//
// Sparse versions of vector operations.  The sparse versions
// detect numerical zeros, so that sparse containers don't store
// numerically zero elements.
//

#if !defined(MPTOOLKIT_LINEARALGEBRA_ITEROPERATIONS_SPARSE_H)
#define MPTOOLKIT_LINEARALGEBRA_ITEROPERATIONS_SPARSE_H

#include "scalar.h"
#include "vectorinterface.h"
#include "vecptriterator.h"
#include <cstring>

template <typename I1, typename I2, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_inner_prod_sparse(I1 i1, I2 i2, Func f, 
                       vector_iterator_ordered, vector_iterator_ordered)
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
   auto NormSq = norm_frob_sq(Result);
   int Count = 1;
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
         result_type Temp = f(*i1, *i2);
         NormSq += norm_frob_sq(Temp);
         ++Count;
         Result += Temp;
         ++i1; ++i2;
      }
   }
   auto Tol = std::numeric_limits<decltype(NormSq)>::epsilon()*10;
   if (norm_frob_sq(Result) < NormSq * Tol * Tol)
      return zero_or_die<value_type>();
   return Result;
}

template <typename I1, typename I2, typename CF, typename Func>
typename make_value_with_zero<typename Func::result_type>::type
iter_coefficient_inner_prod_sparse(I1 i1, I2 i2, CF const& cf, 
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

   // find the first non-zero coefficient
   cf_result x(cf(i2.index()));
   TRACE_IF(norm_frob(x) < 1E-10)(x);
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
      TRACE_IF(norm_frob(x) < 1E-10)(x);
   }

   value_type Result(x * f(*i1, *i2));
   auto NormSq = norm_frob_sq(Result);
   int Count = 1;
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
      {
         result_type Temp = x * f(*i1, *i2);
         NormSq += norm_frob_sq(Temp);
         Result += Temp;
         ++Count;
      }
      ++i1; ++i2;
   }
   auto Tol = std::numeric_limits<decltype(NormSq)>::epsilon()*10;
   if (norm_frob_sq(Result) < NormSq * Tol * Tol)
      return zero_or_die<value_type>();
   return Result;
}

#endif
