// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/coefficient_parallel_prod.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
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
  coefficient_parallel_prod.h

  Created 2005-03-03 Ian McCulloch

  generalized parallel_prod

*/

#if !defined(COEFFICIENT_PARALLEL_PROD_H_HSCKJHUIYH34879YT98374Y)
#define COEFFICIENT_PARALLEL_PROD_H_HSCKJHUIYH34879YT98374Y

#include "matrixoperations.h"

namespace LinearAlgebra
{

// coefficient_parallel_prod

template <typename S, typename T, typename CF,
          typename Nested = Multiplication<typename interface<S>::value_type,
                                           typename interface<T>::value_type>,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct CoefficientParallelProd {};

template <typename S, typename T, typename CF, typename Func>
inline
typename CoefficientParallelProd<S, T, CF, Func>::result_type
coefficient_parallel_prod(S const& x, T const& y, CF const& cf, Func const& f)
{
   return CoefficientParallelProd<S, T, CF, Func>()(x,y,cf,f);
}

template <typename S, typename T, typename CF>
inline
typename CoefficientParallelProd<S, T, CF>::result_type
coefficient_parallel_prod(S const& x, T const& y, CF const& cf)
{
   return CoefficientParallelProd<S, T, CF>()(x,y, cf);
}

template <typename S, typename T, typename CF, typename Func,
          typename Sv, typename Si, typename Tv, typename Ti>
struct CoefficientParallelProd<S, T, CF, Func,
                               VECTOR_EXPRESSION(Sv, Si), VECTOR_EXPRESSION(Tv, Ti)>
   : CoefficientParallelProd<typename EvalExpression<S>::result_type,
                             typename EvalExpression<T>::result_type,
                             CF, Func> {};

template <typename S, typename T, typename CF, typename Func,
          typename Sv, typename Si, typename Tv, typename Ti>
struct CoefficientParallelProd<S, T, CF, Func, LOCAL_VECTOR(Sv, Si), LOCAL_VECTOR(Tv, Ti)>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef CF const& third_argument_type;
   typedef Func const& fourth_argument_type;

   result_type operator()(S const& x, T const& y, CF const& cf) const
   { return iter_coefficient_inner_prod(iterate(x), iterate(y), cf, Func()); }

   result_type operator()(S const& x, T const& y, CF const& cf, Func const& f) const
   { return iter_coefficient_inner_prod(iterate(x), iterate(y), cf, f); }
};

//
// coefficient_parallel_prod_cull
//

template <typename S, typename T, typename CF, typename Float,
          typename Nested = Multiplication<typename interface<S>::value_type,
                                           typename interface<T>::value_type>,
          typename SInterface = typename interface<S>::type,
          typename TInterface = typename interface<T>::type,
          typename Enable = void>
struct CoefficientParallelProdCull {};

template <typename S, typename T, typename CF, typename Float, typename Func>
inline
typename CoefficientParallelProdCull<S, T, CF, Float, Func>::result_type
coefficient_parallel_prod_cull(S const& x, T const& y, CF const& cf, Float const& Tol, Func const& f)
{
   return CoefficientParallelProdCull<S, T, CF, Float, Func>()(x,y,cf,Tol,f);
}

template <typename S, typename T, typename CF, typename Float>
inline
typename CoefficientParallelProdCull<S, T, CF, Float>::result_type
coefficient_parallel_prod_cull(S const& x, T const& y, CF const& cf, Float const& Tol)
{
   return CoefficientParallelProdCull<S, T, CF, Float>()(x,y, cf, Tol);
}

template <typename S, typename T, typename CF, typename Float, typename Func,
          typename Sv, typename Si, typename Tv, typename Ti>
 struct CoefficientParallelProdCull<S, T, CF, Float, Func,
                               VECTOR_EXPRESSION(Sv, Si), VECTOR_EXPRESSION(Tv, Ti)>
   : CoefficientParallelProdCull<typename EvalExpression<S>::result_type,
                             typename EvalExpression<T>::result_type,
                                 CF, Float, Func> {};

template <typename S, typename T, typename CF, typename Float, typename Func,
          typename Sv, typename Si, typename Tv, typename Ti>
struct CoefficientParallelProdCull<S, T, CF, Float, Func, LOCAL_VECTOR(Sv, Si), LOCAL_VECTOR(Tv, Ti)>
{
   typedef typename make_value_with_zero<typename Func::result_type>::type result_type;
   typedef S first_argument_type;
   typedef T second_argument_type;
   typedef CF const& third_argument_type;
   typedef Float const& fourth_argument_type;
   typedef Func const& fifti_argument_type;

   result_type operator()(S const& x, T const& y, CF const& cf, Float const& Tol) const
   { return iter_coefficient_inner_prod_cull(iterate(x), iterate(y), cf, Func(), Tol); }

   result_type operator()(S const& x, T const& y, CF const& cf, Float const& Tol, Func const& f) const
   { return iter_coefficient_inner_prod_cull(iterate(x), iterate(y), cf, f, Tol); }
};

} // namespace LinearAlgebra

#endif
