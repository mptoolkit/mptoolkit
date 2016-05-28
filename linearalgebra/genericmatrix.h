// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/genericmatrix.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
   genericmatrix.h

   Created 2004-05-18 Ian McCulloch

   Top-level interface for matrix classes.  This defines a minimal interface
   for rvalue expressions.  About the only thing you can do with them is assign
   to a concrete matrix class.
*/

#if !defined(GENERICMATRIX_H_HDSJDFHU4YR789Y89ACFH4389H)
#define GENERICMATRIX_H_HDSJDFHU4YR789Y89ACFH4389H

#include "vectorinterface.h"  // for size_type
#include <boost/mpl/bool.hpp>

namespace LinearAlgebra
{

template <typename Derived> 
struct MatrixTraits;

// interface tag for generic matrices
struct GenericMatrixInterface {};

template <typename Scalar, typename Derived>
class GenericMatrix
{
   public:
      typedef Derived                                         derived_type;
      typedef Scalar                                          value_type;
      typedef MatrixTraits<Derived>                           traits_type;
      typedef typename traits_type::expression_type           expression_type;

      typedef GenericMatrixInterface                          interface;
  
      derived_type const& as_derived() const { return static_cast<derived_type const&>(*this); }

      size_type size() const { return this->as_derived().size(); }
      size_type size1() const { return this->as_derived().size1(); }
      size_type size2() const { return this->as_derived().size2(); }

      expression_type eval_expr() const { return this->as_derived().eval_expr(); }

};

template <typename T1, typename D1, typename T2, typename D2>
inline
bool operator==(GenericMatrix<T1, D1> const& x, GenericMatrix<T2, D2> const& y)
{
   return x.eval_expr() == y.eval_expr();
}

template <typename Scalar, typename Derived>
inline
std::ostream& operator<<(std::ostream& out, GenericMatrix<Scalar, Derived> const& M)
{
  return out << M.eval_expr();
}

//
// is_matrix
// 
// traits class to determine if an arbitrary type T is a matrix.
//

namespace Private
{

template <typename T>
T* MakePtr();   // dummy function, not implemented

typedef char IsMatrix[1];
typedef char NotMatrix[2];

NotMatrix& TestIsMatrix(...);

template <typename T, typename Derived>
IsMatrix& TestIsMatrix(GenericMatrix<T, Derived> const*);

template <typename T, int IsMatrix>
struct is_matrix_helper : public boost::mpl::false_
{ };

template <typename T>
struct is_matrix_helper<T, 1> : public boost::mpl::true_
{ };

} // namespace Private   

template <typename T>
struct is_matrix : public Private::is_matrix_helper<T, sizeof(Private::TestIsMatrix(Private::MakePtr<T>()))>
{ };

} // namespace LinearAlgebra

#endif
