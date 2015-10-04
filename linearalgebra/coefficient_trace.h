/* -*- C++ -*- $Id$

  coefficient_trace.h

  Simple implementation of the coefficient_trace function.

  Created 2005-03-18 Ian McCulloch
*/

#if !defined(COEFFICIENT_TRACE_HDSUFHREUHY35879YT893YT879Y)
#define COEFFICIENT_TRACE_HDSUFHREUHY35879YT893YT879Y

#include "matrixoperations.h"

namespace LinearAlgebra
{

template <typename T, typename CF, 
          typename Nested = Trace<typename interface<T>::value_type>, 
          typename Ti = typename interface<T>::type>
struct CoefficientMatrixTrace {};

template <typename T, typename CF>
inline
typename CoefficientMatrixTrace<T, CF>::result_type
coefficient_trace(T const& x, CF const& cf)
{
   return CoefficientMatrixTrace<T, CF>()(x, cf);
}

// implementation

template <typename T, typename CF, typename Nested, typename Tv, typename Ti>
struct CoefficientMatrixTrace<T, CF, Nested, Concepts::SparseMatrix<Tv, Ti>>
{
   typedef typename make_value<typename Nested::result_type>::type result_type;
   typedef T const& first_argument_type;
   typedef CF const& second_argument_type;
   typedef Nested third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type cf,
                          third_argument_type f) const
   {
      typedef typename make_value_with_zero<result_type>::type zval_type;
      zval_type Result = zero<zval_type>();

      typename const_iterator<T>::type I = iterate(x);
      while (I)
      {
         typename const_inner_iterator<T>::type J = iterate(I);
         while (J)
         {
            if (J.index1() == J.index2())
               add(Result, f(*J) * cf(J.index1()));

            ++J;
         }
         ++I;
      }
      return Result;
   }

   result_type operator()(first_argument_type x, second_argument_type cf) const
   {
      return operator()(x, cf, Nested());
   }
};

template <typename T, typename CF, typename Nested, typename Tv, typename Orient, typename Ti>
struct CoefficientMatrixTrace<T, CF, Nested, Concepts::CompressedOuterMatrix<Tv, Orient, Ti>>
{
   typedef typename make_value<typename Nested::result_type>::type result_type;
   typedef T const& first_argument_type;
   typedef CF const& second_argument_type;
   typedef Nested third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type cf,
                          third_argument_type f) const
   {
      typedef typename make_value_with_zero<result_type>::type zval_type;
      zval_type Result = zero<zval_type>();

      for (auto I = iterate(x); I; ++I)
      {
	 auto J = iterate_at(*I, I.index());
	 if (J)
	 {
	    add(Result, f(*J) * cf(J.index()));
	 }
      }
      return Result;
   }

   result_type operator()(first_argument_type x, second_argument_type cf) const
   {
      return operator()(x, cf, Nested());
   }
};

template <typename T, typename CF, typename Nested, typename Tv, typename Orient, typename Ti>
struct CoefficientMatrixTrace<T, CF, Nested, Concepts::DenseMatrix<Tv, Orient, Ti>>
{
   typedef typename make_value<typename Nested::result_type>::type result_type;
   typedef T const& first_argument_type;
   typedef CF const& second_argument_type;
   typedef Nested third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type cf,
                          third_argument_type f) const
   {
      typedef typename make_value_with_zero<result_type>::type zval_type;
      zval_type Result = zero<zval_type>();

      for (auto I = iterate(x); I; ++I)
      {
	 auto J = iterate_at(*I, I.index());
	 if (J)
	 {
	    add(Result, f(*J) * cf(J.index1()));
	 }
      }
      return Result;
   }

   result_type operator()(first_argument_type x, second_argument_type cf) const
   {
      return operator()(x, cf, Nested());
   }
};

template <typename T, typename CF, typename Nested, typename Tv, typename Ti>
struct CoefficientMatrixTrace<T, CF, Nested, Concepts::DiagonalMatrix<Tv, Ti>>
{
   typedef typename make_value<typename Nested::result_type>::type result_type;
   typedef T const& first_argument_type;
   typedef CF const& second_argument_type;
   typedef Nested third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type cf,
                          third_argument_type f) const
   {
      typedef typename make_value_with_zero<result_type>::type zval_type;
      zval_type Result = zero<zval_type>();

      for (auto I = iterate(x.diagonal()); I; ++I)
      {
	 add(Result, f(*I) * cf(I.index()));
      }
      return Result;
   }

   result_type operator()(first_argument_type x, second_argument_type cf) const
   {
      return operator()(x, cf, Nested());
   }
};

} // namespace LinearAlgebra

#endif
