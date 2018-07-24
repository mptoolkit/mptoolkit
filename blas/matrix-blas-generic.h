// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-blas-generic.icc
//
// Copyright (C) 2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// generic versions of BLAS functions, for BLAS-like extensions and
// blas functions that work for non-scalar types

#if !defined(MPTOOLKIT_BLAS_MATRIX_BLAS_GENERIC_H)
#define MPTOOLKIT_BLAS_MATRIX_BLAS_GENERIC_H

#include "functors.h"

namespace blas
{

//
// vector
//

template <typename T>
void vector_fill(T const& alpha, int N, T* y, int incy);

template <typename T, typename U>
void
vector_copy(int M, T const* x, int incx, U* y, int incy);

template <typename T, typename U, typename V>
void
vector_copy_scaled(int N, T const& alpha, U const* A, int lda, V* B, int ldb);

template <typename  T>
void vector_conj(int N, T* x, int incx);

template <typename T, typename U>
void
vector_add(int N, T const* x, int incx, U* y, int incy);

template <typename T, typename U, typename V>
void
vector_add_scaled(int N, T const& alpha, U const* x, int incx, V* y, int incy);

template <typename T, typename U>
void
vector_scale(int N, T const& alpha, U* y, int incy);

template <typename T>
void
vector_sum(int N, T const* x, int incx, T& r);

#if 0
template <typename T>
void
vector_inner_prod(int N, double const* x, int incx, double const* y, int incy, double& r);

void
vector_inner_prod(int N, std::complex<double> const* x, int incx, 
		  std::complex<double> const* y, int incy,
		  std::complex<double>& r);
#endif

template <typename T, typename U, typename V, typename Nested>
inline
void
vector_inner_prod_nested(int N, T const* x, int incx, U const* y, int incy, V& r, Nested&& n);

template <typename T, typename U, typename V, typename Nested>
void
vector_add_inner_prod_nested(int N, T const* x, int incx, T const* y, int incy, V& r, Nested&& n);

//
// matrix
//

//void matrix_clear(int M, int N, double* A, int lda);
//void matrix_clear(int N, std::complex<double>* y, int incy);

template <typename T>
void matrix_fill(T const& alpha, int M, int N, T* A, int lda);

template <typename T>
void matrix_conj(int M, int N, T* A, int lda);

template <typename T>
decltype(norm_frob_sq(std::declval<T>()))
matrix_norm_frob_sq(int M, int N, T const* A, int lda);

template <typename T, typename U>
void
matrix_copy(int M, int N, T const* A, int lda, U* B, int ldb);

template <typename T, typename U>
void
matrix_copy(char Atrans, int M, int N, T const* A, int lda, U* B, int ldb);

template <typename T, typename U, typename V>
void
matrix_copy_scaled(char Atrans, int M, int N, T alpha, U const* A, int lda, V* B, int ldb);

template <typename T, typename U>
void
matrix_add(char Atrans, int M, int N, T const* A, int lda, U* , int ldbB);

template <typename T, typename U, typename V>
void
matrix_add_scaled(char Atrans, int M, int N, T alpha, U const* A, int lda, V* B, int ldb);

#if 0
template <typename T, typename U, typename V>
void
matrix_inner_prod(char Atrans, char Btrans, int M, int N,
		  double const* A, int ldA,
		  double const* B, int ldB,
		  double& r);

void
matrix_inner_prod(char Atrans, char Btrans, int M, int N,
		  std::complex<double> const* A, int ldA,
		  std::complex<double> const* B, int ldB,
		  std::complex<double>& r);
#endif

template <typename T, typename U, typename V, typename Nested>
void
matrix_inner_prod_nested(char Atrans, char Btrans, int M, int N,
			 T const* A, int ldA,
			 U const* B, int ldB,
			 V& r, Nested&& n);

// level 3

// We constrain beta to be a scalar here
template <typename T, typename U, typename V, typename W, typename X>
void
gemm(char Atrans, char Btrans, int M, int N, int K, T const& alpha,
     U const* A, int lda,
     V const* B, int ldb,
     W beta, X* C, int ldc);

// product diagonal with general matrix
template <typename Xt, typename Bt, typename Ct>
void
dgmm(int M, int K, 
     Xt const* x, int incx,
     Bt const* B, int ldb,
     Ct* C, int ldc);

template <typename Xt, typename Ct>
void
dgmm_inplace(int M, int K,
	     Xt const* x, int incx,
	     Ct* B, int ldc);

// product general matrix with diaqonal matrix
template <typename At, typename Yt, typename Ct>
void
gdmm(int M, int N,
     At const* A, int lda,
     Yt const* y, int incy,
     Ct* C, int ldc);

template <typename At, typename Yt>
void
gdmm_inplace(int M, int N,
	     At* A, int lda,
	     Yt const* y, int incy);

} // namespace blas

#include "matrix-blas-generic.icc"

#endif
