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

#if !defined(MPTOOLKIT_BLAS_MATRIX_LOWLEVEL_REFERENCE_H)
#define MPTOOLKIT_BLAS_MATRIX_LOWLEVEL_REFERENCE_H

#include "matrix-lowlevel-base.h"

namespace blas
{

inline char const* matrix_blas_library()
{
   return "Standard BLAS (no extensions in use)";
};

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

void
matrix_copy(char Atrans, int M, int N, double const* A, int lda, double* , int ldbB);

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

} // namespace blas

#include "matrix-lowlevel-reference.icc"

#endif

