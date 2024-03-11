// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/coefficient_operations.h
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
  coefficient_operations.h

  Created 2005-03-03 Ian McCulloch

  Generalizations of transpose and matrix multiplication that
  introduce a scale factor.

  Rough first attempt, no expression templates.
*/

#if !defined(COEFFICIENT_OPERATIONS_H_HSCKJHUIYH34879YT98374Y)
#define COEFFICIENT_OPERATIONS_H_HSCKJHUIYH34879YT98374Y

#include "matrixoperations.h"
#include "coefficient_multiply.h"
#include "coefficient_multiply_cull.h"
#include <boost/utility/result_of.hpp>

namespace LinearAlgebra
{

//
// CoefficientMatrixProductProxy
//

template <typename T1, typename T2, typename CoeffFunc, typename NestedMult>
class CoefficientMatrixProductProxy;

template <typename T1, typename T2, typename CoeffFunc, typename NestedMult>
struct abstract_interface<CoefficientMatrixProductProxy<T1, T2, CoeffFunc, NestedMult> >
{
   typedef typename matrix_abstract_or<abstract_interface<T1>,
                                       abstract_interface<T2> >::type type;
};

template <typename T1, typename T2, typename CoeffFunc, typename NestedMult>
class CoefficientMatrixProductProxy
{
   public:
      BOOST_MPL_ASSERT((is_const_proxy_reference<T1>));
      BOOST_MPL_ASSERT((is_const_proxy_reference<T2>));

      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      typedef T1 reference1;
      typedef T2 reference2;

      typedef typename basic_type<T1>::type matrix1_type;
      typedef typename basic_type<T2>::type matrix2_type;

      typedef CoeffFunc coefficient_functor_type;

      typedef NestedMult functor_type;

      size_type size1() const { return Size1<matrix1_type>()(x_); }
      size_type size2() const { return Size2<matrix2_type>()(y_); }

      reference1 matrix1() const { return x_; }
      reference2 matrix2() const { return y_; }

      CoefficientMatrixProductProxy(reference1 x, reference2 y, CoeffFunc const& cf,
                         functor_type const& f = functor_type()) : x_(x), y_(y), cf_(cf), f_(f)
      {
         using LinearAlgebra::size1; using LinearAlgebra::size2;
         DEBUG_PRECONDITION_EQUAL(size2(x), size1(y));
      }

      functor_type functor() const { return f_; }
      coefficient_functor_type coefficient_functor() const { return cf_; }

   private:
      CoefficientMatrixProductProxy& operator=(CoefficientMatrixProductProxy const&);
   // not implemented

      reference1 x_;
      reference2 y_;

      coefficient_functor_type cf_;
      functor_type f_;
};

// interface

template <typename T1, typename T2, typename CoeffFunc, typename Nested>
struct interface<CoefficientMatrixProductProxy<T1, T2, CoeffFunc, Nested> >
{
   typedef typename interface<T1>::value_type value1_type;
   typedef typename interface<T2>::value_type value2_type;
   //   typedef typename make_value<typename Nested::result_type>::type value_type;
   typedef typename make_value<
      typename boost::result_of<
         Nested(value1_type, value2_type)
      >::type
   >::type value_type;

   using type = Concepts::MatrixExpression<value_type,
                                           CoefficientMatrixProductProxy<T1, T2, CoeffFunc, Nested>>;
};

// Assign

template <typename LHS, typename T1, typename T2, typename CF, typename F>
struct AssignExpression<LHS&, CoefficientMatrixProductProxy<T1, T2, CF, F> >
{
   typedef LHS& result_type;
   typedef LHS& first_argument_type;
   typedef CoefficientMatrixProductProxy<T1, T2, CF, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      assign_coefficient_product2(x, y.matrix1(), y.matrix2(),
                                  y.coefficient_functor(), y.functor());
      return x;
   }
};

// Add

template <typename LHS, typename T1, typename T2, typename CF, typename F>
struct AddExpression<LHS&, CoefficientMatrixProductProxy<T1, T2, CF, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef CoefficientMatrixProductProxy<T1, T2, CF, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      add_coefficient_product2(x, y.matrix1(), y.matrix2(),
                               y.coefficient_functor(), y.functor());
   }
};

// Subtract

template <typename LHS, typename T1, typename T2, typename CF, typename F>
struct SubtractExpression<LHS&, CoefficientMatrixProductProxy<T1, T2, CF, F> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef CoefficientMatrixProductProxy<T1, T2, CF, F> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      subtract_coefficient_product2(x, y.matrix1(), y.matrix2(),
                                    y.coefficient_functor(), y.functor());
   }
};

//
// CoefficientMatrixMatrixMultiplication
//

template <typename M1, typename M2,
          typename CF,
          typename Nested = Multiplication<typename interface<M1>::value_type,
                                           typename interface<M2>::value_type>,
          typename M1i = typename interface<M1>::type,
          typename M2i = typename interface<M2>::type>
struct CoefficientMatrixMatrixMultiplication {};

template <typename M1, typename M2, typename CF>
inline
typename CoefficientMatrixMatrixMultiplication<M1, M2, CF>::result_type
coefficient_multiply(M1 const& m1, M2 const& m2, CF const& cf)
{
   return CoefficientMatrixMatrixMultiplication<M1, M2, CF>()(m1, m2, cf);
}

template <typename M1, typename M2, typename CF, typename F>
inline
typename CoefficientMatrixMatrixMultiplication<M1, M2, CF, F>::result_type
coefficient_multiply(M1 const& m1, M2 const& m2, CF const& cf, F const& f)
{
   return CoefficientMatrixMatrixMultiplication<M1, M2, CF, F>()(m1, m2, cf, f);
}

template <typename M1, typename M2,
          typename CF,
          typename Nested,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct CoefficientMatrixMatrixMultiplication<M1, M2, CF, Nested,
                                             Concepts::MatrixExpression<M1v, M1i>,
                                             Concepts::MatrixExpression<M2v, M2i>>
{
   typedef M1 const& first_argument_type;
   typedef M2 const& second_argument_type;
   typedef CF const& third_argument_type;
   typedef Nested const& fourth_argument_type;
   typedef CoefficientMatrixProductProxy<typename make_const_reference<M1>::type,
                                         typename make_const_reference<M2>::type,
                                         CF, Nested> result_type;

   result_type operator()(M1 const& m1, M2 const& m2,
                          CF const& cf) const { return result_type(m1, m2, cf); }

   result_type operator()(M1 const& m1, M2 const& m2,
                          CF const& cf, Nested const& f) const
      { return result_type(m1, m2, cf, f); }
};

//
// coefficient_transpose
//

template <typename MatType, typename Func>
typename make_value<typename Transpose<MatType>::result_type>::type
coefficient_transpose(MatType const& M, Func Coefficient)
{
   typedef typename make_value<typename Transpose<MatType>::result_type>::type
      result_type;

   result_type Result = transpose(M);

   typename iterator<result_type>::type I = iterate(Result);
   while (I)
   {
      typename inner_iterator<result_type>::type J = iterate(I);
      while (J)
      {
         *J *= Coefficient(J.index2(), J.index1());
         ++J;
      }
      ++I;
   }
   return Result;
}

template <typename MatType, typename Func, typename Nest>
typename make_value<typename MatrixTranspose<MatType, Nest>::result_type>::type
coefficient_transpose(MatType const& M, Func Coefficient, Nest n)
{
   typedef typename make_value<typename MatrixTranspose<MatType, Nest>::result_type>::type
      result_type;

   result_type Result = transpose(M, n);

   typename iterator<result_type>::type I = iterate(Result);
   while (I)
   {
      typename inner_iterator<result_type>::type J = iterate(I);
      while (J)
      {
         *J *= Coefficient(J.index2(), J.index1());
         ++J;
      }
      ++I;
   }
   return Result;
}

// prod_cull()

template <typename T1, typename T2, typename CoeffFunc, typename NestedMult, typename Float>
class CoefficientMatrixProductCullProxy;

template <typename T1, typename T2, typename CoeffFunc, typename NestedMult, typename Float>
struct abstract_interface<CoefficientMatrixProductCullProxy<T1, T2, CoeffFunc, NestedMult, Float> >
{
   typedef typename matrix_abstract_or<abstract_interface<T1>,
                                       abstract_interface<T2> >::type type;
};

template <typename T1, typename T2, typename CoeffFunc, typename NestedMult, typename Float>
class CoefficientMatrixProductCullProxy
{
   public:
      BOOST_MPL_ASSERT((is_const_proxy_reference<T1>));
      BOOST_MPL_ASSERT((is_const_proxy_reference<T2>));

      // mark this type as a const proxy reference
      typedef boost::mpl::true_ const_proxy;

      typedef T1 reference1;
      typedef T2 reference2;

      typedef typename basic_type<T1>::type matrix1_type;
      typedef typename basic_type<T2>::type matrix2_type;

      typedef CoeffFunc coefficient_functor_type;

      typedef NestedMult functor_type;

      size_type size1() const { return Size1<matrix1_type>()(x_); }
      size_type size2() const { return Size2<matrix2_type>()(y_); }

      reference1 matrix1() const { return x_; }
      reference2 matrix2() const { return y_; }

      CoefficientMatrixProductCullProxy(reference1 x, reference2 y, CoeffFunc const& cf,
                                        Float const& Tol,
                                        functor_type const& f = functor_type())
         : x_(x), y_(y), cf_(cf), Tol_(Tol), f_(f)
      {
         using LinearAlgebra::size1; using LinearAlgebra::size2;
         DEBUG_PRECONDITION_EQUAL(size2(x), size1(y));
      }

      functor_type functor() const { return f_; }
      coefficient_functor_type coefficient_functor() const { return cf_; }
      Float const& tol() const { return Tol_; }

   private:
      CoefficientMatrixProductCullProxy& operator=(CoefficientMatrixProductCullProxy const&);
   // not implemented

      reference1 x_;
      reference2 y_;

      coefficient_functor_type cf_;
      Float Tol_;
      functor_type f_;
};

// interface

template <typename T1, typename T2, typename CoeffFunc, typename Nested, typename Float>
struct interface<CoefficientMatrixProductCullProxy<T1, T2, CoeffFunc, Nested, Float> >
{
   typedef typename interface<T1>::value_type value1_type;
   typedef typename interface<T2>::value_type value2_type;
   //   typedef typename make_value<typename Nested::result_type>::type value_type;
   typedef typename make_value<
      typename boost::result_of<
         Nested(value1_type, value2_type)
      >::type
   >::type value_type;

   using type = Concepts::MatrixExpression<value_type,
                                           CoefficientMatrixProductCullProxy<T1, T2,
                                                                             CoeffFunc, Nested, Float>>;
};

// Assign

template <typename LHS, typename T1, typename T2, typename CF, typename F, typename Float>
struct AssignExpression<LHS&, CoefficientMatrixProductCullProxy<T1, T2, CF, F, Float> >
{
   typedef LHS& result_type;
   typedef LHS& first_argument_type;
   typedef CoefficientMatrixProductCullProxy<T1, T2, CF, F, Float> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      try_resize(x, size1(y), size2(y));
      assign_coefficient_cull_product2(x, y.matrix1(), y.matrix2(),
                                       y.coefficient_functor(), y.functor(), y.tol());
      return x;
   }
};

// Add

template <typename LHS, typename T1, typename T2, typename CF, typename F, typename Float>
struct AddExpression<LHS&, CoefficientMatrixProductCullProxy<T1, T2, CF, F, Float> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef CoefficientMatrixProductCullProxy<T1, T2, CF, F, Float> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      add_coefficient_cull_product2(x, y.matrix1(), y.matrix2(),
                                    y.coefficient_functor(), y.functor(), y.tol());
   }
};

// Subtract

template <typename LHS, typename T1, typename T2, typename CF, typename F, typename Float>
struct SubtractExpression<LHS&, CoefficientMatrixProductCullProxy<T1, T2, CF, F, Float> >
{
   typedef void result_type;
   typedef LHS& first_argument_type;
   typedef CoefficientMatrixProductCullProxy<T1, T2, CF, F, Float> const& second_argument_type;
   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      PRECONDITION_EQUAL(size1(x), size1(y));
      PRECONDITION_EQUAL(size2(x), size2(y));
      subtract_coefficient_product2(x, y.matrix1(), y.matrix2(),
                                    y.coefficient_functor(), y.functor(), y.tol());
   }
};

//
// CoefficientMatrixMatrixMultiplicationCull
//

template <typename M1, typename M2,
          typename CF, typename Float,
          typename Nested = Multiplication<typename interface<M1>::value_type,
                                           typename interface<M2>::value_type>,
          typename M1i = typename interface<M1>::type,
          typename M2i = typename interface<M2>::type>
struct CoefficientMatrixMatrixMultiplicationCull {};

template <typename M1, typename M2, typename CF, typename Float>
inline
typename CoefficientMatrixMatrixMultiplicationCull<M1, M2, CF, Float>::result_type
coefficient_multiply_cull(M1 const& m1, M2 const& m2, CF const& cf, Float const& Tol)
{
   return CoefficientMatrixMatrixMultiplicationCull<M1, M2, CF, Float>()(m1, m2, cf, Tol);
}

template <typename M1, typename M2, typename CF, typename Float, typename F>
inline
typename CoefficientMatrixMatrixMultiplicationCull<M1, M2, CF, Float, F>::result_type
coefficient_multiply_cull(M1 const& m1, M2 const& m2, CF const& cf, Float const& Tol, F const& f)
{
   return CoefficientMatrixMatrixMultiplicationCull<M1, M2, CF, Float, F>()(m1, m2, cf, Tol, f);
}

template <typename M1, typename M2,
          typename CF, typename Float,
          typename Nested,
          typename M1v, typename M1i,
          typename M2v, typename M2i>
struct CoefficientMatrixMatrixMultiplicationCull<M1, M2, CF, Float, Nested,
                                                 Concepts::MatrixExpression<M1v, M1i>,
                                                 Concepts::MatrixExpression<M2v, M2i>>
{
   typedef M1 const& first_argument_type;
   typedef M2 const& second_argument_type;
   typedef CF const& third_argument_type;
   typedef Float const& fourth_argument_type;
   typedef Nested const& fifth_argument_type;
   typedef CoefficientMatrixProductCullProxy<typename make_const_reference<M1>::type,
                                             typename make_const_reference<M2>::type,
                                             CF, Nested, Float> result_type;

   result_type operator()(M1 const& m1, M2 const& m2,
                          CF const& cf, Float const& Tol) const { return result_type(m1, m2, cf, Tol); }

   result_type operator()(M1 const& m1, M2 const& m2,
                          CF const& cf, Float const& Tol, Nested const& f) const
   { return result_type(m1, m2, cf, Tol, f); }
};

} // namespace LinearAlgebra

#endif
