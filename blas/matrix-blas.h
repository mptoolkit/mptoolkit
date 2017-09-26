// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-blas.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_MATRIX_BLAS_H)
#define MPTOOLKIT_BLAS_MATRIX_BLAS_H

#include "matrix.h"
#include "matrix-lowlevel-blas.h"

namespace blas
{

template <typename T, typename U, typename V>
inline
void gemm(T alpha, BlasMatrix<T, Matrix<T>, U> const& A, BlasMatrix<T, Matrix<T>, V> const& B,
          T beta, Matrix<T>& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   gemm(A.trans(), B.trans(), A.rows(), A.cols(), B.cols(), alpha, data(A.as_derived()),
	A.leading_dimension(), data(B.as_derived()), B.leading_dimension(), beta,
        C.data(), C.leading_dimension());
}

template <typename T, typename U>
inline
void matrix_copy_scaled(T alpha, BlasMatrix<T, Matrix<T>, U> const& A, Matrix<T> const& C)
{
   matrix_copy_scaled(A.trans(), A.rows(), A.cols(), alpha, A.leading_dimension(), data(A.as_derived()),
                      C.data(), C.leading_dimension());
}

template <typename T, typename U>
inline
void matrix_copy(BlasMatrix<T, Matrix<T>, U> const& A, Matrix<T> const& C)
{
   matrix_copy(A.trans(), A.rows(), A.cols(), A.leading_dimension(), data(A.as_derived()),
                      C.data(), C.leading_dimension());
}

} // namespace blas

#endif
