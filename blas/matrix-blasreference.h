// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-lowlevel-reference.h
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

#if !defined(MPTOOLKIT_BLAS_MATRIX_BLASREFERENCE_H)
#define MPTOOLKIT_BLAS_MATRIX_BLASREFERENCE_H

#include "matrix-blas-base.h"

namespace blas
{

inline char const* matrix_blas_library()
{
   return "Standard BLAS (no extensions in use)";
};

//
// vector
//

void
vector_copy(int M, double const* x, int incx, double* y, int incy);

void
vector_copy(int N, std::complex<double> const* x, int incx, std::complex<double>* y, int incy);

void
vector_copy_scaled(int N, double alpha, double const* A, int lda, double* B, int ldb);

void
vector_copy_scaled(int N, std::complex<double> alpha,
                   std::complex<double> const* x, int incx,
                   std::complex<double>* y, int incy);

void
vector_add(int N, double const* x, int incx, double* y, int incy);

void
vector_add(int N, std::complex<double> const* x, int incx,
           std::complex<double>* y, int incy);

void
vector_add_scaled(int N, double alpha, double const* x, int incx, double* y, int incy);

void
vector_add_scaled(int N, std::complex<double> alpha,
                  std::complex<double> const* x, int incx,
                  std::complex<double>* y, int incy);

//
// matrix
//

void
matrix_copy(char Atrans, int M, int N, double const* x, int incx, double* y, int incy);

void
matrix_copy(char Atrans, int M, int N,
            std::complex<double> const* A, int lda,
            std::complex<double>* B, int ldb);

void
matrix_copy_scaled(char Atrans, int M, int N, double alpha, double const* A, int lda, double* B, int ldb);

void
matrix_copy_scaled(char Atrans, int M, int N, std::complex<double> alpha,
                   std::complex<double> const* A, int lda,
                   std::complex<double>* B, int ldb);

void
matrix_add(char Atrans, int M, int N, double const* A, int lda, double* , int ldbB);

void
matrix_add(char Atrans, int M, int N,
           std::complex<double> const* A, int lda,
           std::complex<double>* B, int ldb);

void
matrix_add_scaled(char Atrans, int M, int N, double alpha, double const* A, int lda, double* B, int ldb);

void
matrix_add_scaled(char Atrans, int M, int N, std::complex<double> alpha,
                  std::complex<double> const* A, int lda,
                  std::complex<double>* B, int ldb);

// level 3

inline
void
gemm(char Atrans, char Btrans, int M, int N, int K, std::complex<double> alpha,
     std::complex<double> const* A, int lda,
     std::complex<double> const* B, int ldb,
     std::complex<double> beta, std::complex<double>* C, int ldc)
{
   CHECK(Atrans != 'R')("R trans is not yet implemented!");
   CHECK(Btrans != 'R')("R trans is not yet implemented!");
   BLAS::zgemm(Atrans, Btrans, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
}

} // namespace blas

#include "matrix-blasreference.icc"

#endif

