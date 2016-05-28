// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/coefficient_parallel_prod.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  coefficient_parallel_prof.h

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

} // namespace LinearAlgebra

#endif
