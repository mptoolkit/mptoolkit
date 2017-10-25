// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-eigen.h
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

//
// A simple dense matrix class designed for scalar types.
//

#if !defined(MPTOOLKIT_BLAS_MATRIX_EIGEN_H)
#define MPTOOLKIT_BLAS_MATRIX_EIGEN_H

#include "blas/vector.h"
#include "blas/matrix.h"

namespace blas
{

//
// DiagonalizeSymmetric
//
// Diagonalizes a matrix in-place.  The matrix is replaced by the transform matrix, with
// the eigenvectors as sucessive column-vectors.  For input matrix M,
// X = M' is the transform matrix, E is the diagonal matrix of eigenvalues,
// we have MX = XE
//

inline
Vector<double> DiagonalizeSymmetric(Matrix<double>& M);

template <typename U>
void DiagonalizeSymmetric(Matrix<double>& M, BlasVector<double, U, cpu_tag>& v)

// Version that takes a proxy-reference
template <typename U>
void DiagonalizeSymmetric(Matrix<double>& M, BlasVector<double, U, cpu_tag>&& v)

inline
Vector<double>
DiagonalizeHermitian(Matrix<double>& M)
{
   return DiagonalizeSymmetric(M);
}

template <typename U>
inline
void DiagonalizeHermitian(Matrix<double>& M, BlasVector<double, U, cpu_tag>& v)
{
   DiagonalizeSymmeric(M, v.as_derived());
}

template <typename U>
inline
void DiagonalizeHermitian(Matrix<double>& M, BlasVector<double, U, cpu_tag>&& v)
{
   DiagonalizeSymmeric(M, static_cast<U&&>(v.as_derived()));
}

//
// complex
//

inline
Vector<double> DiagonalizeHermitian(Matrix<std::complex<double>>& M);

template <typename U>
void DiagonalizeHermitian(Matrix<std::complex<double>>& M, BlasVector<double, U, cpu_tag>& v);

template <typename U>
void DiagonalizeHermitian(Matrix<std::complex<double>>& M, BlasVector<double, U, cpu_tag>&& v);

} // namespace blas

#include "matrix-eigen.icc"

#endif
