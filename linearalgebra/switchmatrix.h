// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/switchmatrix.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
/* -*- C++ -*- $Id: scalarmatrix.h 1149 2012-04-18 03:12:37Z ianmcc $

  switchmatrix.h

  A SwitchMatrix is a matrix that is a variant of one of
  ScalarMatrix<complex<double> >
  Matrix<double>
  Matrix<complex<double> >

  The usual operations are defined, that attempt as much as possible
  to keep the most efficient representation, ScalarMatrix first,
  then Matrix<double>, and lastly Matrix<complex<double> >

  Created 2014-04-17 Ian McCulloch
*/

#if !defined(MPTOOLKIT_LINEARALGEBRA_SWITCHMATRIX_H)
#define MPTOOLKIT_LINEARALGEBRA_SWITCHMATRIX_H

#include "matrix.h"
#include "scalarmatrix.h"
#include <boost/variant.hpp>

namespace LinearAlgebra
{

typedef ScalarMatrix<std::complex<double> > ComplexScalarMatrix;
typedef Matrix<double>                      DoubleMatrix;
typedef Matrix<st::complex<double> >        ComplexMatrix;

typedef boost::variant<ComplexScalarMatrix, DoubleMatrix, ComplexMatrix> SwitchVariantType;

class SwitchMatrix
{
   public:
      typedef std::complex<double> data_type;

      typedef boost::mpl::false_ proxy;
      typedef boost::mpl::false_ const_proxy;
      typedef boost::mpl::false_ immediate;

      typedef data_type value_type;

      SwitchMatrix();

      SwitchMatrix(SwitchVariantType const& x) : Data(x) {}

      SwitchMatrix(size_type s1, size_type s2);

      SwitchMatrix(size_type s1, size_type s2, value_type x);

      size_type size1() const;
      size_type size2() const;

      void resize(int r, int c);

      void fill(value_type x);

      SwitchMatrix& operator=(SwitchMatrix const& m);

      template <typename U>
      typename boost::enable_if<is_matrix<U>, SwtichMatrix&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_matrix<U>, SwtichMatrix&>::type
      operator=(NoAliasProxy<U> const& x);

      template <typename RHS>
      SwitchMatrix& operator+=(RHS const& m)
      {
         return add_copy(*this, m);
      }

      template <typename RHS>
      SwitchMatrix& operator-=(RHS const& m)
      {
         return subtract_copy(*this, m);
      }

      SwitchVariantType const& data() const { return Data; }
      SwitchVariantType& data() { return Data; }

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v) const
      {
         return boost::apply_visitor(v, Data);
      }

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v)
      {
         return boost::apply_visitor(v, Data);
      }

   private:
      SwitchVariantType Data;
};

template <typename  Func>
inline
typename Func::return_type
apply_visitor(Func const& f, SwitchMatrix const& m)
{
   return boost::apply_visitor(f, m.data());
}

template <typename  Func>
inline
typename Func::return_type
apply_visitor(Func const& f, SwitchMatrix& m)
{
   return boost::apply_visitor(f, m.data());
}

// interface

struct interface<SwitchMatrix>
{
   typedef SwitchMatrix::value_type value_type;
   typedef RECTANGULAR_MATRIX(value_type, SwitchMatrix) type;
};

// iterators are not defined for SwitchMatrix

// GetMatrixElement is not defined for SwitchMatrix

// Conj

struct SwitchVariantTypeConj : public boost::static_visitor<SwitchVariantType>
{
   ComplexScalarMatrix operator()(ComplexScalarMatrix const& x) const
   {
      return ComplexScalarMatrix(x.size1(), conj(x.value()));
   }

   DoubleMatrix const& operator()(DoubleMatrix const& x) const
   {
      return x;
   }

   ComplexMatrix operator()(ComplexMatrix const& x) const
   {
      return conj(x);
   }
};

struct Conj<SwitchMatrix>
{
   typedef SwitchMatrix const& argument_type;
   typedef SwitchMatrix result_type;

   result_type operator()(argument_type x) const
   {
      return SwitchMatrix(apply_visitor(SwitchVariantTypeConj(), x));
   }
};

// Herm

struct SwitchVariantTypeHerm : public boost::static_visitor<SwitchVariantType>
{
   ComplexScalarMatrix operator()(ComplexScalarMatrix const& x) const
   {
      return ComplexScalarMatrix(x.size1(), conj(x.value()));
   }

   DoubleMatrix const& operator()(DoubleMatrix const& x) const
   {
      return herm(x);
   }

   ComplexMatrix operator()(ComplexMatrix const& x) const
   {
      return herm(x);
   }
};

struct Herm<SwitchMatrix>
{
   typedef SwitchMatrix const& argument_type;
   typedef SwitchMatrix result_type;

   result_type operator()(argument_type x) const
   {
      return SwitchMatrix(apply_visitor(SwitchVariantTypeHerm(), x));
   }
};

// Multiply

// This is the operator*= for a SwitchMatrix with another SwitchMatrix.
// The return type is a SwitchVariantType, but we are allowed to modify the
// left hand operator too

template <typename U>
struct SwitchVariantTypeMultiply;

struct SwitchVariantTypeMultiply<double> : public boost::static_visitor<SwitchVariantType>
{
   template <typename T>
   T& operator()(T& x, double y)
   {
      x *= y;
      return x;
   }
}

struct SwitchVariantTypeMultiply<std::complex<double> >
   : public boost::static_visitor<SwitchVariantType>
{
   ScalarMatrix& operator()(ScalarMatrix& x, std::complex<double> y)
   {
      x.value() *= y;
      return x;
   }

   SwitchVariantType operator()(DoubleMatrix& x, std::complex<double> y)
   {
      if (y.imag() == 0)
      {
         x *= y.real();
         return x;
      }
      // else
      return ComplexMatrix(x * y.value());
   }

   ComplexMatrix& operator()(ComplexMatrix& x, std::complex<double> y)
   {
      x *= y;
      return x;
   }
}

template <>
struct Multiply<SwitchMatrix&, double>
{
   typedef SwitchMatrix& result_type;
   typedef SwitchMatrix& first_argument_type;
   typedef U second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      return apply_visitor(SwitchVariantTypeMultiply<double>(y), x);
   }
};

template <>
struct Multiply<SwitchMatrix&, std::complex<double> >
{
   typedef SwitchMatrix& result_type;
   typedef SwitchMatrix& first_argument_type;
   typedef U second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      return apply_visitor(SwitchVariantTypeMultiply<std::complex<double> >(y), x);
   }
};

struct SwitchVariantTypeMultiplySelf : public boost::static_visitor<SwitchVariantType>
{
   ScalarMatrix& operator()(ScalarMatrix& x, ScalarMatrix const& y)
   {
      DEBUG_CHECK_EQUAL(x.size2(), y.size1());
      x.value() *= y.value();
      return x;
   }

   SwitchVariantType operator()(DoubleMatrix& x, ScalarMatrix const& y)
   {
      DEBUG_CHECK_EQUAL(x.size2(), y.size1());
      if (y.value().imag() == 0)
      {
         x *= y.value().real();
         return x;
      }
      // else
      return ComplexMatrix(x * y.value());
   }

   ComplexMatrix& operator()(ComplexMatrix& x, ScalarMatrix const& y)
   {
      DEBUG_CHECK_EQUAL(x.size2(), y.size1());
      x *= y.value();
      return x;
   }

   SwitchVariantType operator()(ScalarMatrix& x, DoubleMatrix const& y)
   {
      DEBUG_CHECK_EQUAL(x.size2(), y.size1());
      if (y.value()imag() == 0)
      {
         return DoubleMatrix(x.value().real() * y);
      }
      // else
      return ComplexMatrix(x.value() * y);
   }

   DoubleMatrix& operator()(DoubleMatrix& x, DoubleMatrix const& y)
   {
      x *= y;
      return x;
   }

   ComplexMatrix& operator()(ComplexMatrix& x, DoubleMatrix const& y)
   {
      x *= y;
      return x;
   }

   ComplexMatrix& operator()(ScalarMatrix& x, ComplexMatrix const& y)
   {
      DEBUG_CHECK_EQUAL(x.size2(), y.size1());
      return x.value() * y;
   }

   ComplexMatrix& operator()(DoubleMatrix& x, ComplexMatrix const& y)
   {
      return x * y;
   }

   ComplexMatrix& operator()(ComplexMatrix& x, ComplexMatrix const& y)
   {
      x *= y;
      return x;
   }
};


template <typename U>
struct Multiply<SwitchMatrix&, SwitchMatrix>
{
   typedef SwitchMatrix& result_type;
   typedef SwitchMatrix& first_argument_type;
   typedef U second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      x = apply_visitor(SwitchVariantTypeMultiplySelf(), x, y);
      return x;
   }
};

// we don't need to further specalize Multiply, since it defaults to calling
// assign_copy(x, x*y) which is fine here.


// Transform

template  <typename F>
struct SwitchVariantTypeTransform : public boost::static_visitor<SwitchVariantType>
{
   SwitchVariantTypeTransform(F const& f_) : f(f_) {}
   template <typename T>
   SwitchVariantType operator()(T const& m) const
   {
      return transform(m, f);
   }
};

template <typename F>
struct TransformMatrix<SwitchMatrix, F>
{
   typedef SwitchMatrix const& argument_type;
   typedef argument_type first_argument_type;
   typedef F const& second_argument_type;

   typedef SwitchMatrix result_type;

   result_type operator()(argument_type e) const
   {
      return apply_visitor(SwitchVariantTypeTransform<F>(), e);
   }

   result_type operator()(argument_type e, second_argument_type f) const
   {
      return apply_visitor(SwitchVariantTypeTransform<F>(f), e);
   }
};

// non-const Transform isn't defined


// BinaryTransform isn't defined

// multiplication

template <typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<SwitchMatrix, SwitchMatrix, F,
                                  RECTANGULAR_MATRIX(Sv, Si), RECTANGULAR_MATRIX(Tv, Ti)>
{
   typedef typename is_commutative<F>::type commutative;

   typedef typename F::result_type FwdResult;
   typedef typename make_value<FwdResult>::type result_value_type;

   typedef SwitchMatrix<result_value_type> result_type;
   typedef SwitchMatrix<S> const& first_argument_type;
   typedef SwitchMatrix<T> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return result_type(x.size1(), F()(x.value(), y.value()), cdirect());
   }
};

template <typename T, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<SwitchMatrix, T, F,
                                  RECTANGULAR_MATRIX(Sv, Si), ANY_MATRIX(Tv, Ti)>
{
   typedef SwitchMatrixMultiplication<Sv, T, F> Fwd;
   typedef typename Fwd::result_type result_type;
   typedef SwitchMatrix<S> const& first_argument_type;
   typedef T const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return Fwd()(x.value(), y);
   }
};

template <typename S, typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixMatrixMultiplication<S, SwitchMatrix, F,
                                  ANY_MATRIX(Sv, Si), RECTANGULAR_MATRIX(Tv, Ti)>
{
   typedef MatrixScalarMultiplication<S, Tv, F> Fwd;
   typedef typename Fwd::result_type result_type;
   typedef S const& first_argument_type;
   typedef SwitchMatrix<T> const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y)
   {
      DEBUG_PRECONDITION_EQUAL(x.size2(), y.size1());
      return Fwd()(x, y.value());
   }
};

//
// DirectProduct
//

template <typename F, typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixDirectProduct<SwitchMatrix, SwitchMatrix, F,
                           DIAGONAL_MATRIX(Sv, Si), DIAGONAL_MATRIX(Tv, Ti)>
{
   typedef typename F::result_type ValType;
   typedef typename make_value<ValType>::type result_value_type;
   typedef SwitchMatrix<result_value_type> result_type;

   typedef SwitchMatrix<S> const& first_argument_type;
   typedef SwitchMatrix<T> const& second_argument_type;
   typedef F const& third_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return result_type(x.size1() * y.size1(), F()(x.value(), y.value()), cdirect());
   }

   result_type operator()(first_argument_type x, second_argument_type y,
                          third_argument_type f) const
   {
      return result_type(x.size1() * y.size1(), f(x.value(), y.value()), cdirect());
   }
};

// InnerProd

template <>
struct InnerProd<SwitchMatrix, SwitchMatrix>
{
   typedef SwitchMatrix const& first_argument_type;
   typedef SwitchMatrix const& second_argument_type;
   typedef std::complex<double> result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.size1(), y.size1());
      return inner_prod(x.value(), y.value()) * result_type(x.size1());
   }
};

template <typename S,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<S, SwitchMatrix, InnerProd<Sv, Tv>, ANY_MATRIX(Sv, Si),
                       DIAGONAL_MATRIX(Tv, Ti)>
{
   typedef S const& first_argument_type;
   typedef SwitchMatrix const& second_argument_type;
   typedef typename result_value<InnerProd<Sv, Tv> >::type result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return inner_prod(trace(x), y.value());
   }
};

template <typename T,
          typename Sv, typename Si, typename Tv, typename Ti>
struct MatrixInnerProd<SwitchMatrix, T, InnerProd<Sv, Tv>, DIAGONAL_MATRIX(Sv, Si),
                       ANY_MATRIX(Tv, Ti)>
{
   typedef SwitchMatrix const& first_argument_type;
   typedef T const& second_argument_type;
   typedef typename result_value<InnerProd<Sv, Tv> >::type result_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      return inner_prod(x.value(), trace(y));
   }
};



} // namespace LinearAlgebra

#endif
