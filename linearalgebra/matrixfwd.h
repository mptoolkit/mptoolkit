// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixfwd.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  forward declarations for some matrix types that are used as temporaries.

  Created 2005-03-10 Ian McCulloch
*/

#if !defined(MATRIXFWD_H_FHURH34Y89Y98Y9OY89PQYW9)
#define MATRIXFWD_H_FHURH34Y89Y98Y9OY89PQYW9

#include "matrixinterface.h"

namespace LinearAlgebra
{

template <typename Scalar, typename Orientation = RowMajor>
class Matrix;

template <typename T, typename Orientation = RowMajor,
          typename InnerType = MapVector<T>, typename OuterType = Vector<InnerType> >
class SparseMatrix;

template <typename T, typename DiagonalType = Vector<T> >
class DiagonalMatrix;

template <typename Scalar>
class ScalarMatrix;

template <typename T, typename I = typename interface<T>::type>
struct make_matrix_from_interface {};

// forward make_value_from_interface for matrix types to
// make_matrix_from_interface
template <typename T, typename Sv, typename Si>
 struct make_value_from_interface<T, Concepts::AnyMatrix<Sv, Si>>
   : make_matrix_from_interface<T, Concepts::AnyMatrix<Sv, Si>> {};

template <typename M, typename T, typename Ti>
struct make_matrix_from_interface<M, Concepts::LocalMatrix<T, Ti>>
{
   typedef Matrix<T> type;
};

template <typename M, typename T, typename Orient, typename Ti>
struct make_matrix_from_interface<M, Concepts::DenseMatrix<T, Orient, Ti>>
{
   typedef Matrix<T, Orient> type;
};

template <typename M, typename T, typename Ti>
struct make_matrix_from_interface<M, Concepts::SparseMatrix<T, Ti>>
{
   typedef SparseMatrix<T> type;
};

template <typename M, typename T, typename Orient, typename Ti>
struct make_matrix_from_interface<M, Concepts::CompressedOuterMatrix<T, Orient, Ti>>
{
   typedef SparseMatrix<T, Orient> type;
};

template <typename M, typename T, typename Ti>
struct make_matrix_from_interface<M, Concepts::DiagonalMatrix<T, Ti>>
{
   typedef DiagonalMatrix<T> type;
};

template <typename M, typename T, typename Ti>
struct make_matrix_from_interface<M, Concepts::ScalarMatrix<T, Ti>>
{
   typedef ScalarMatrix<T> type;
};

//
// Expresion types
// Here, we forward make_matrix_from_interface to make_matrix_from_abstract,
// and we use the abstract_interface.

template <typename T, typename AI = typename abstract_interface<T>::type>
struct make_matrix_from_abstract {};

template <typename T>
struct make_matrix_from_abstract<T, matrix_abstract_dense>
{
   using type = Matrix<T>;
};

template <typename T>
struct make_matrix_from_abstract<T, matrix_abstract_sparse>
{
   using type = SparseMatrix<T>;
};

template <typename T, typename Tv, typename Ti>
struct make_matrix_from_interface<T, Concepts::MatrixExpression<Tv, Ti>>
   : make_matrix_from_abstract<Tv, typename abstract_interface<T>::type> {};



} // namespace LinearAlgebra

#endif
