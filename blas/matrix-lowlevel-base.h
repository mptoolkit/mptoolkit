// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/blas-base.h
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

// standard BLAS functions

#if !defined(MPTOOLKIT_BLAS_MATRIX_LOWLEVEL_BASE_H)
#define MPTOOLKIT_BLAS_MATRIX_LOWLEVEL_BASE_H

#include "common/blas1f.h"
#include "common/blas2f.h"
#include "common/blas3f.h"

namespace blas
{

// gemm functions

inline
void
gemm(char Atrans, char Btrans, int M, int N, int K, double alpha, double const* A, int lda,
     double const* B, int ldb, double beta, double* C, int ldc)
{
   BLAS::dgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline
void
vector_copy(int N, double const* x, int incx, double* y, int incy)
{
   BLAS::dcopy(N, x, incx, y, incy);
}

inline
void
vector_copy(int N, std::complex<double> const* x, int incx, std::complex<double>* y, int incy)
{
   BLAS::zcopy(N, x, incx, y, incy);
}

} // namespace blas

#endif
