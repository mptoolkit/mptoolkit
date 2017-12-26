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

template <typename I1, typename I2, typename T>
inline
void vector_copy_from_stl(I1 start, I2 finish, T* y, int incy)
{
   while (start != finish)
   {
      *y = *start;
      y += incy;
      ++start;
   }
}

void vector_clear(int N, double* y, int incy);

void vector_clear(int N, std::complex<double>* y, int incy);

void vector_fill(double alpha, int N, double* y, int incy);

void vector_fill(std::complex<double> alpha, int N, std::complex<double>* y, int incy);

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

void
vector_sum(int N, double const* x, int incx, double& r);

void
vector_inner_prod(int N, double const* x, int incx, double const* y, int incy, double& r);

void
vector_inner_prod(int N, std::complex<double> const* x, int incx, 
		  std::complex<double> const* y, int incy,
		  std::complex<double>& r);

template <typename Nested>
inline
void
vector_inner_prod_nested(int N, double const* x, int incx, double const* y, int incy, double& r, Nested&&)
{
   vector_inner_prod(N, x, incx, y, incy, r);
}

template <typename Nested>
void
vector_inner_prod_nested(int N, std::complex<double> const* x, int incx, 
			 std::complex<double> const* y, int incy,
			 std::complex<double>& r, Nested&&)
{
   vector_inner_prod(N, x, incx, y, incy, r);
}

template <typename T>
std::enable_if_t<blas::is_numeric_v<T>, void>
vector_add_inner_prod(int N, T const* x, int incx, T const* y, int incy, T& r);

template <typename T, typename Nested>
std::enable_if_t<blas::is_numeric_v<T>, void>
vector_add_inner_prod_nested(int N, T const* x, int incx, T const* y, int incy, T& r, Nested&&)
{
   vector_add_inner_prod(N, x, incx, y, incy, r);
}

//
// matrix
//

void matrix_clear(int M, int N, double* A, int lda);
void matrix_clear(int N, std::complex<double>* y, int incy);

void matrix_fill(double alpha, int M, int N, double* A, int lda);
void matrix_fill(std::complex<double> alpha, int M, int N, std::complex<double>* A, int lda);

template <typename T>
decltype(norm_frob_sq(std::declval<T>()))
matrix_norm_frob_sq(int M, int N, T const* A, int lda);

void
matrix_copy(char Atrans, int M, int N, double const* A, int lda, double* B, int ldb);

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

#if 0
void
matrix_add_inner_prod(char Atrans, char Btrans, int M, int N,
		      double const* A, int ldA,
		      double const* B, int ldB,
		      double& r);

void
matrix_add_inner_prod(char Atrans, char Btrans, int M, int N,
		      std::complex<double> const* A, int ldA,
		      std::complex<double> const* B, int ldB,
		      std::complex<double>& r);
#endif

template <typename T, typename Nested>
std::enable_if_t<std::is_floating_point<T>::value, void>
matrix_inner_prod_nested(char Atrans, char Btrans, int M, int N,
			 T const* A, int ldA,
			 T const* B, int ldB,
			 T& r, Nested&)
{
   matrix_inner_prod(Atrans, Btrans, M, N, A, ldA, B, ldB, r);
}

template <typename T, typename Nested>
std::enable_if_t<std::is_floating_point<T>::value, void>
matrix_inner_prod_nested(char Atrans, char Btrans, int M, int N,
			 std::complex<T> const* A, int ldA,
			 std::complex<T> const* B, int ldB,
			 std::complex<T>& r, Nested&)
{
   matrix_inner_prod(Atrans, Btrans, M, N, A, ldA, B, ldB, r);
}

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

