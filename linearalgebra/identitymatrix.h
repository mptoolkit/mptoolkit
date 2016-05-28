// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/identitymatrix.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(IDENTITYMATRIX_H_CJDHCJKHUIYT389Y987Y34897GYO8)
#define IDENTITYMATRIX_H_CJDHCJKHUIYT389Y987Y34897GYO8

#include "matrix.h"

namespace LinearAlgebra
{

template <typename T>
   class IdentityMatrix : public GenericMatrix<T, IdentityMatrix<T> >
{
   public:
      typedef T value_type;

      IdentityMatrix();

      explicit IdentityMatrix(size_type s1);
      IdentityMatrix(size_type s1, size_type s2);

      size_type size1() const { return size_; }
      size_type size2() const { return size_; }

      IdentityMatrix eval_expr() const { return *this; }

   private:
      size_type size_;
};

template <typename T>
struct MatrixTraits<IdentityMatrix<T> >
{
   typedef IdentityMatrix<T> expression_type;
};

template <typename S, typename T>
struct MatrixMatrixProduct<IdentityMatrix<S>, IdentityMatrix<T> >
{
   typedef typename BinaryOperator<Multiplication, S, T>::value_type ValType;
   typedef IdentityMatrix<ValType> result_type;
   typedef IdentityMatrix<ValType> value_type;

   typedef IdentityMatrix<S> first_argument_type;
   typedef IdentityMatrix<T> second_argument_type;

   result_type operator()(first_argument_type const& x, second_argument_type const& y)
   {
      DEBUG_PRECONDITION(x.size2() == y.size1());
      return result_type(x.size1(), y.size2());
   }
};

template <typename S>
struct Conj<IdentityMatrix<S> >
{
   typedef IdentityMatrix<S> argument_type;
   typedef IdentityMatrix<S> result_type;
   typedef result_type value_type;

   result_type operator()(argument_type const& x) const { return x; }
};

#if 0
template <typename S>
struct is_identity_function<Conj<IdentityMatrix<S> > > : public boost::mpl::true_ { };
#endif

} // namespace LinearAlgebra

#endif
