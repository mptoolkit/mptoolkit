// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/blas2f.h
//
// Copyright (C) 2001-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  blas2f.h

  C++ interface to BLAS level 2

  Created 2001-04-09 Ian McCulloch

  Ideas from LAPACK++
*/

#if !defined(MPTOOLKIT_COMMON_BLAS2F_H)
#define MPTOOLKIT_COMMON_BLAS2F_H

#include "fortran.h"
#include "restrict.h"

namespace BLAS
{

using namespace Fortran;

void sgemv(char trans, integer M, integer N, float alpha, 
	   float const* A, integer lda, float const* dx, 
	   integer incx, float beta, float* restrict dy, integer incy);


void dgemv(char trans, integer M, integer N, double alpha, 
	   double const* A, integer lda, double const* dx, 
	   integer incx, double beta, double* restrict dy, integer incy);

void dgbmv(char trans, integer M, integer N, integer kl, 
	   integer ku, double alpha, double const* A, integer lda, 
	   double const* dx, integer incx, double beta, 
	   double* restrict dy, integer incy);

void dsymv(char uplo, integer N, double alpha, double const* A, 
	   integer lda, double const* dx, integer incx, double beta,
	   double* restrict dy, integer incy);

void dsbmv(char uplo, integer N, integer k, double alpha, 
	   double const* A, integer lda, double const* dx, 
	   integer incx, double beta, double* restrict dy, integer incy);

void dspmv(char uplo, integer N, double alpha, double* restrict AP, 
	   double const* dx, integer incx, double beta, double* restrict dy, 
	   integer incy);

void dtrmv(char uplo, char trans, char diag, integer N, 
	   double const* A, integer lda, double* restrict dx, 
	   integer incx);

void dtrsv(char uplo, char trans, char diag, integer N, 
	   double const* A, integer lda, double* restrict dx, integer incx);

void dger(integer M, integer N, double alpha, 
	  double const* dx, integer incx, double const* dy, integer incy, 
	  double* restrict A, integer lda);

void dsyr(char uplo, integer N, double alpha, double const* dx, 
	  integer incx, double* restrict A, integer lda);

void dspr(char uplo, integer N, double alpha, double const* dx, 
	  integer incx, double* restrict AP);

void dsyr2(char uplo, integer N, double alpha, double* const dx, 
	   integer incx, double* const dy, integer incy, double* restrict A, 
	   integer lda);

void dspr2(char uplo, integer N, double alpha, double const* dx, 
	   integer incx, double const* dy, integer incy, double* restrict AP);


namespace raw
{

extern "C"
{
     void F77NAME(sgemv)(char* trans, integer* M, integer* N, float* alpha, 
			 const float* A, integer* lda, const float* dx, 
			 integer* incx, float* beta, float* dy, integer* incy);


     void F77NAME(dgemv)(char* trans, integer* M, integer* N, double* alpha, 
			 const double* A, integer* lda, const double* dx, 
			 integer* incx, double* beta, double* dy, integer* incy);

     void F77NAME(dgbmv)(char* trans, integer* M, integer* N, integer* kl, 
			 integer* ku, double* alpha, const double* A, integer* lda, 
			 const double* dx, integer* incx, double* beta, 
			 double* dy, integer* incy);

     void F77NAME(dsymv)(char* uplo, integer* N, double* alpha, const double* A, 
			 integer* lda, const double* dx, integer* incx, double* beta,
			 double* dy, integer* incy);

     void F77NAME(dsbmv)(char* uplo, integer* N, integer* k, double* alpha, 
			 const double* A, integer* lda, const double* dx, 
			 integer* incx, double* beta, double* dy, integer* incy);

     void F77NAME(dspmv)(char* uplo, integer* N, double* alpha, double* AP, 
			 double const* dx, integer const* incx, double* beta, double* dy, 
			 integer* incy);

     void F77NAME(dtrmv)(char* uplo, char* trans, char* diag, const integer* N, 
			 const double* A, integer* lda, const double* dx, 
			 integer* incx);

     void F77NAME(dtrsv)(char* uplo, char* trans, char* diag, const integer* N, 
			 double* A, integer* lda, double* dx, integer* incx);

     void F77NAME(dger)(integer* M, integer* N, double* alpha, 
			double* dx, integer* incx, double* dy, integer* incy, 
			double* A, integer* lda);

     void F77NAME(dsyr)(char* uplo, integer* N, double* alpha, double* dx, 
			integer* incx, double* A, integer* lda);

     void F77NAME(dspr)(char* uplo, integer* N, double* alpha, double* dx, 
			integer* incx, double* AP);

     void F77NAME(dsyr2)(char* uplo, integer* N, double* alpha, double* dx, 
			 integer* incx, double* dy, integer* incy, double* A, 
			 integer* lda);

     void F77NAME(dspr2)(char* uplo, integer* N, double* alpha, double* dx, 
			 integer* incx, double* dy, integer* incy, double* AP);

} // extern "C"

} // namespace raw

// inlines

// float

inline void sgemv(char trans, integer M, integer N, float alpha, 
		  float const* A, integer lda, float const* dx, 
		  integer incx, float beta, float* restrict dy, integer incy)
{
   raw::F77NAME(sgemv)(&trans, &M, &N, &alpha, A, &lda, dx, &incx, &beta, dy, &incy);
}

// double

inline void dgemv(char trans, integer M, integer N, double alpha, 
		  double const* A, integer lda, double const* dx, 
		  integer incx, double beta, double* restrict dy, integer incy)
{
   raw::F77NAME(dgemv)(&trans, &M, &N, &alpha, A, &lda, dx, &incx, &beta, dy, &incy);
}

inline void dgbmv(char trans, integer M, integer N, integer kl, 
		  integer ku, double alpha, double const* A, integer lda, 
		  double const* dx, integer incx, double beta, 
		  double* restrict dy, integer incy)
{
   raw::F77NAME(dgbmv)(&trans, &M, &N, &kl, &ku, &alpha, A, &lda, dx, &incx, &beta, dy, &incy);
}

inline void dsymv(char uplo, integer N, double alpha, double const* A, 
		  integer lda, double const* dx, integer incx, double beta,
		  double* restrict dy, integer incy)
{
   raw::F77NAME(dsymv)(&uplo, &N, &alpha, A, &lda, dx, &incx, &beta, dy, &incy);
}

inline void dsbmv(char uplo, integer N, integer k, double alpha, 
		  double const* A, integer lda, double const* dx, 
		  integer incx, double beta, double* restrict dy, integer incy)
{
   raw::F77NAME(dsbmv)(&uplo, &N, &k, &alpha, A, &lda, dx, &incx, &beta, dy, &incy);
}

#if 0
inline void dspmv(char uplo, integer N, double alpha, double* restrict AP, 
		  double const* dx, integer incx, double beta, double* restrict dy, 
		  integer incy)
{
   raw::F77NAME(dspmv)(&uplo, &N, &alpha, AP, dx, &incx, &beta, dy, &incy);
}

inline void dtrmv(char uplo, char trans, char diag, integer N, 
		  double const* A, integer lda, double* restrict dx, 
		  integer incx)
{
   raw::F77NAME(dtrmv)(&uplo, &trans, &diag, &N, A, &lda, dx, &incx);
}

inline void dtrsv(char uplo, char trans, char diag, integer N, 
		  double const* A, integer lda, double* restrict dx, integer incx)
{
   raw::F77NAME(dtrsv)(&uplo, &trans, &diag, &N, A, &lda, dx, &incx);
}

inline void dger(integer M, integer N, double alpha, 
		 double const* dx, integer incx, double const* dy, integer incy, 
		 double* restrict A, integer lda)
{
   raw::F77NAME(dger)(&M, &N, &alpha, dx, &incx, &dy, &incy, A, &lda);
}

inline void dsyr(char uplo, integer N, double alpha, double const* dx, 
		 integer incx, double* restrict A, integer lda)
{
   raw::F77NAME(dsyr)(&uplo, &N, &alpha, dx, &incx, A, &lda);
}

inline void dspr(char uplo, integer N, double alpha, double const* dx, 
		 integer incx, double* restrict AP)
{
   raw::F77NAME(dspr)(&uplo, &N, *alpha, dx, &incx, AP);
}

inline void dsyr2(char uplo, integer N, double alpha, double* const dx, 
		  integer incx, double* const dy, integer incy, double* restrict A, 
		  integer lda)
{
   raw::F77NAME(dsyr2)(&uplo, &N, &alpha, dx, &incx, dy, &incy, A, &lda);
}

inline void dspr2(char uplo, integer N, double alpha, double const* dx, 
		  integer incx, double const* dy, integer incy, double* restrict AP)
{
   raw::F77NAME(dspr2)(&uplo, &N, &alpha, dx, &incx, dy, &incy, AP);
}
#endif

} // namespace BLAS

/*
BLAS documentation

**************************************************************************************


sgemv, dgemv, cgemv, zgemv 
Matrix-vector product for a general matrix

FORMAT  
  {S,D,C,Z}GEMV (trans, m, n, alpha, a, lda, x, incx, beta, y, incy)

Arguments  

  trans               character*1
                      On entry, specifies the operation to be performed:

                      If trans = 'N' or 'n', the operation is y  =  alpha*Ax
                      + beta*y.

                      If trans = 'T' or 't', the operation is y  =
                      alpha*transp(A)*x + beta*y.

                      If trans = 'C' or 'c', the operation is y  =
                      alpha*conjug_transp(A)*x + beta*y.
                      On exit, trans is unchanged.

  m                   integer*4
                      On entry, the number of rows of the matrix A; m >= 0.
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns of the matrix A; n >=
                      0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n. The leading m by n part of the array contains the
                      elements of the matrix A.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      MAX(1,m).
                      On exit, lda is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array containing the vector
                      x.  When trans is equal to 'N' or (1+(n-1)*|incx|).
                      Otherwise, the length is at least (1+(m-1)*|incx|).
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar beta.
                      On exit, beta is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array containing the vector
                      x.  When trans is equal to 'N' or (1+(m-1)*|incy|).
                      Otherwise, the length is at least (1+(n-1)*|incy|).

                      If beta= 0, y need not be set.  If beta is not equal to
                      zero, the incremented array Y must contain the vector
                      y.
                      On exit, y is overwritten by the updated vector y.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

Description  
  The _GEMV subprograms compute a matrix-vector product for either a general
  matrix or its transpose: y  =  alpha*Ax + beta*y
    y  =  alpha*transp(A)*x + beta*y

  In addition to these operations, the CGEMV and ZGEMV subprograms compute
  the matrix-vector product for the conjugate transpose:
    y  =  alpha*conjug_transp(A)*x + beta*y

  alpha and beta are scalars, x and y are vectors, and A is an m by n matrix.

Example
  
  REAL*8 A(20,20), X(20), Y(20), alpha, beta
  INCX = 1
  INCY = 1
  LDA = 20
  M = 20
  N = 20
  alpha = 1.0D0
  beta = 0.0D0
  CALL DGEMV('T',M,N,alpha,A,LDA,X,INCX,beta,Y,INCY)

  This FORTRAN code computes the product y = transp(A)*x.

  COMPLEX*8 A(20,20), X(20), Y(20), alpha, beta
  INCX = 1
  INCY = 1
  LDA = 20
  M = 20
  N = 20
  alpha = (1.0, 1.0)
  beta = (0.0, 0.0)
  CALL CGEMV('T',M,N,alpha,A,LDA,X,INCX,beta,Y,INCY)

  This FORTRAN code computes the product y = transp(A)*x.


**************************************************************************************


sgbmv, ddbmv, cgbmv, zgbmv 
Matrix-vector product for a general band matrix

FORMAT  
  {S,D,C,Z}GBMV (trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)

Arguments  

  trans               character*1
                      On entry, specifies the operation to be performed:

                      If trans = 'N' or 'n', the operation is y  =  alpha*Ax
                      + beta*y.

                      If trans = 'T' or 't', the operation is y  =
                      alpha*transp(A)*x + beta*y.

                      If trans = 'C' or 'c', the operation is y  =
                      alpha*conjug_transp(A)*x + beta*y.
                      On exit, trans is unchanged.

  m                   integer*4
                      On entry, the number of rows of the matrix A; m >= 0.
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns of the matrix A; n >=
                      0.
                      On exit, n is unchanged.

  kl                  integer*4
                      On entry, the number of sub-diagonals of the matrix A;
                      kl >= 0.
                      On exit, kl is unchanged.

  ku                  integer*4
                      On entry, the number of super-diagonals of the matrix
                      A; ku >= 0.
                      On exit, ku is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n. The leading m by n part of the array contains the
                      elements of the matrix A, supplied column by column.
                      The leading diagonal of the matrix is stored in row (ku
                      + 1) of the array, the first super-diagonal is stored
                      in row ku starting at position 2,  the first sub-
                      diagonal is stored in row (ku + 2) starting at position
                      1,  and so on. Elements in the array A that do not
                      correspond to elements in the matrix (such as the top
                      left ku by ku triangle) are not referenced.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      (kl+ku+1).
                      On exit, lda is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array containing the vector
                      x.  When trans is equal to 'N' or (1+(n-1)*|incx|).
                      Otherwise, the length is at least (1+(m-1)*|incx|).
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar beta.
                      On exit, beta is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array containing the vector
                      x.  When trans is equal to 'N' or (1+(m-1)*|incy|).
                      Otherwise, the length is at least (1+(n-1)*|incy|).

                      If beta= 0, y need not be set. If beta is not equal to
                      zero, the incremented array Y must contain the vector
                      y.
                      On exit, y is overwritten by the updated vector y.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

Description  
  The _GBMV subprograms compute a matrix-vector product for either a general
  band matrix or its transpose: y  =  alpha*Ax + beta*y
    y  =  alpha*transp(A)*x + beta*y

  In addition to these operations, the CGBMV and ZGBMV subprograms compute a
  matrix-vector product for the conjugate transpose:
    y  =  alpha*conjug_transp(A)*x + beta*y

  alphaand betaare scalars, x and y are vectors, and A is an m by n band
  matrix.

Example 
 
  COMPLEX*16 A(5,20), X(20), Y(20), alpha, beta
  M = 5
  N = 20
  KL = 2
  KU = 2
  alpha = (1.0D0, 2.0D0)
  LDA = 5
  INCX = 1
  beta = (0.0D0, 0.0D0)
  INCY = 1
  CALL ZGBMV('N',M,N,KL,KU,alpha,A,LDA,X,INCX,beta,Y,INCY)

  This FORTRAN code multiplies a pentadiagonal matrix A by the vector x to
  get the vector y.  The operation is y  =  Ax.  where A is stored in banded
  storage form.


**************************************************************************************


ssymv, dsymv, chemv, zhemv 
Matrix-vector product for a symmetric or hermitian matrix

FORMAT  
  {S,D}SYMV (uplo, n, alpha, a, lda, x, incx, beta, y, incy) {C,Z}HEMV (uplo,
  n, alpha, a, lda, x, incx, beta, y, incy)

Arguments  
  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the array A is referenced:

                      If uplo = 'U' or 'u', the upper-triangular part of A is
                      referenced.

                      If uplo = 'L' or 'l', the lower-triangular part of A is
                      referenced.
                      On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n.

  When uplo specifies the upper portion of the matrix, the leading n by n
  part of the array contains the upper-triangular part of the matrix, and the
  lower-triangular part of array A is not referenced.

  When uplo specifies the lower  portion of the matrix,  the leading n by n
  part of the array contains the lower-triangular part of the matrix, and the
  upper-triangular part of array A is not referenced.

  For CHEMV and ZHEMV routines,  the imaginary parts of the diagonal elements
  are not accessed, need not be  set, and are assumed to be zero.
  On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      MAX(1,n).
                      On exit, lda is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar beta.
                      On exit, beta is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).

  If beta= 0, y need not be set.  If betais not equal to zero, the
  incremented array Y must contain the vector y.
  On exit, y is overwritten by the updated vector y.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

Description  
  SSYMV and DSYMV compute a matrix-vector product for a real symmetric
  matrix.  CHEMV and ZHEMV compute a matrix-vector product for a complex
  Hermitian matrix.  Both products are described by the following operation:
  y  = alpha*Ax + beta*y

  alpha and beta are scalars, x and y are vectors with n elements, and A is
  an n by n matrix. In the case of SSYMV and DSYMV, matrix A is a symmetric
  matrix and in the case of CHEMV and ZHEMV, matrix A is a Hermitian matrix.

Example  

  REAL*8 A(100,40), X(40), Y(40), alpha, beta
  N = 40
  INCX = 1
  INCY = 1
  alpha = 1.0D0
  beta = 0.0D0
  LDA = 100
  CALL DSYMV('U',N,alpha,A,LDA,X,INCX,beta,Y,INCY)

  This FORTRAN code computes the product y  =  Ax where A is a symmetric
  matrix, of order 40, with its upper-triangular part stored.

  COMPLEX*8 A(100,40), X(40), Y(40), alpha, beta
  N = 40
  INCX = 1
  INCY = 1
  alpha = (1.0, 0.5)
  beta = (0.0, 0.0)
  LDA = 100
  CALL CHEMV('U',N,alpha,A,LDA,X,INCX,beta,Y,INCY)

  This FORTRAN code computes the product y  =  Ax where A is a Hermitian
  matrix, of order 40, with its upper-triangular part stored.


**************************************************************************************


ssbmv, dsbmv, chbmv, zhbmv 
Matrix-vector product for a symmetric or hermitian band matrix

FORMAT  
  {S,D}SBMV (uplo, n, k, alpha, a, lda, x, incx, beta, y, incy) {C,Z}HBMV
  (uplo, n, k, alpha, a, lda, x, incx, beta, y, incy)

Arguments
  
  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the array A is referenced:

     If uplo = 'U' or 'u', the upper-triangular part of A is referenced.

     If uplo = 'L' or 'l', the lower-triangular part of A is referenced.
     On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  k                   integer*4
                      On entry, if uplo specifies the upper portion of matrix
                      A, k represents the number of super-diagonals of the
                      matrix. If uplo specifies the lower portion, k is the
                      number of subdiagonals; k >= 0.
                      On exit, k is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n.

  When uplo specifies the upper portion of the matrix, the leading (k + 1) by
  n part of the array must contain the upper-triangular band part of the
  matrix, supplied column by column. The main diagonal of the matrix is
  stored in row (k + 1) of the array, the first super-diagonal is stored in
  row k starting at position 2, and so on.  The top left k by k triangle of
  the array A is not referenced.

  When uplo specifies the lower portion of the matrix, the leading (k + 1) by
  n part of the array must contain the lower-triangular band part of the
  matrix, supplied column by column. The main diagonal of the matrix is
  stored in row 1 of the array, the first sub-diagonal is stored in row 2,
  starting at position 1, and so on. The bottom right k by k triangle of the
  array A is not referenced.

  For CHBMV and ZHBMV routines,  the imaginary parts of the diagonal elements
  are not accessed, need not be  set, and are assumed to be zero.
  On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >= (k+1).
                      On exit, lda is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar beta.
                      On exit, beta is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).

  If beta= 0, y need not be set.  If betais not equal to zero, the
  incremented array Y must contain the vector y.
  On exit, y is overwritten by the updated vector y.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

Description  
  SSBMV and DSBMV compute a matrix-vector product for a real symmetric band
  matrix. CHBMV and ZHBMV compute a matrix-vector product for a complex
  Hermitian band matrix. Both products are described by the following
  operation: y  = alpha*Ax + beta*y

  alphaand betaare scalars, and x and y are vectors with n elements. In the
  case of SSBMV and DSBMV, A is a symmetric matrix and in the case of CHBMV
  and ZHBMV, A is a Hermitian matrix.

Example  

  REAL*8 A(2,10), X(10), Y(10), alpha, beta
  N = 10
  K = 1
  alpha = 2.0D0
  LDA = 2
  INCX = 1
  beta = 1.0D0
  INCY = 1
  CALL DSBMV('U',N,K,alpha,A,LDA,X,INCX,beta,Y,INCY)

  This FORTRAN code computes the product y  =  alpha*Ax + y) where A is a
  symmetric tridiagonal matrix, with A stored in upper-triangular form.

  COMPLEX*8 A(2,10), X(10), Y(10), alpha, beta
  N = 10
  K = 1
  alpha = (2.0, 2.2)
  LDA = 2
  INCX = 1
  beta = (1.0, 0.0)

  This FORTRAN code computes the product y  =  alpha*Ax + y) where A is a
  Hermitian tridiagonal matrix, with the upper diagonal of A stored.


**************************************************************************************


sspmv, dspmv, chpmv, zhpmv 
Matrix-vector product for a symmetric or hermitian matrix stored in packed form

FORMAT  
  {S,D}SPMV (uplo, n, alpha, ap, x, incx, beta, y, incy) {C,Z}HPMV (uplo, n,
  alpha, ap, x, incx, beta, y, incy)

Arguments  
  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the matrix A is supplied in the
                      packed array AP:

                      If uplo = 'U' or 'u', the upper-triangular part of A is
                      supplied.

                      If uplo = 'L' or 'l', the lower-triangular part of A is
                      supplied.
                      On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  ap                  real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array AP of length at least
                      n(n + 1)/2.

                      If uplo specifies the upper triangular part of the
                      matrix A, the array contains those elements of the
                      matrix, packed sequentially, column by column, so that
                      AP(1) contains a(11), AP(2) and AP(3) contain a(12) and
                      a(22) respectively, and so on.

                      If uplo specifies the lower triangular part to the
                      matrix A, the array contains those elements of the
                      matrix, also packed sequentially, so that AP(1)
                      contains a(11), AP(2) and AP(3) contain a(21) and a(31)
                      respectively, and so on.

  For CHPMV and ZHPMV routines, the imaginary parts of the diagonal elements
  are not accessed, need not be set, and are assumed to be zero.
  On exit, ap is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar beta.
                      On exit, beta is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).

                      If beta= 0, y need not be set.  If beta is not equal to
                      zero, the incremented array Y must contain the vector
                      y.
                      On exit, y is overwritten by the updated vector y.

  incy                integer*4
                      On entry, the increment  for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

Description  
  SSPMV and DSPMV compute a matrix-vector product for a real symmetric matrix
  stored in packed form. CHPMV and ZHPMV compute a matrix-vector product for
  a complex Hermitian matrix stored in packed form. Both products are
  described by the following operation: y  = alpha*Ax + beta*y

  alpha and beta are scalars, and x and y are vectors with n elements. A is
  an n by n matrix. In the case of SSPMV and DSPMV, matrix A is a symmetric
  matrix and in the case of CHPMV and ZHPMV, matrix A is a Hermitian matrix.

Example  

  COMPLEX*16 AP(250), X(20), Y(20), alpha, beta
  N = 20
  alpha = (2.3D0, 8.4D0)
  INCX = 1
  beta = (4.0D0, 3.3D0)
  INCY = 1
  CALL ZHPMV('L',N,alpha,AP,X,INCX,beta,Y,INCY)

  This FORTRAN code computes the product y  =  alpha*Ax + beta*y where A is a
  Hermitian matrix with its lower-triangular part stored in packed form in
  AP.


**************************************************************************************


strmv, dtrmv, ctrmv, ztrmv 
Marix-vector product for a triangular matrix

FORMAT  
  {S,D,C,Z}TRMV (uplo, trans, diag, n, a, lda, x, incx)

Arguments  
  uplo                character*1
                      On entry, specifies whether the matrix A is an upper-
                      or lower-triangular matrix:

                      If uplo = 'U' or 'u', A is an upper-triangular matrix.

                      If uplo = 'L' or
                      On exit, uplo is unchanged.

  trans               character*1
                      On entry, specifies the operation to be performed:

                      If trans = 'N' or 'n', the operation is y  =  alpha*Ax
                      + beta*y.

                      If trans = 'T' or 't', the operation is y  =
                      alpha*transp(A)*x + beta*y.

                      If trans = 'C' or 'c', the operation is y  =
                      alpha*conjug_transp(A)*x + beta*y.
                      On exit, trans is unchanged.
                      
  diag                character*1
                      On entry, specifies whether the matrix A is unit-
                      triangular:

                      If diag = 'U' or 'u', A is a unit-triangular matrix.

                      If diag = 'N' or 'n', A is not a unit-triangular
                      matrix.
                      On exit, diag is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n.

  When uplo specifies the upper portion of the matrix, the leading n by n
  part of the array contains the upper-triangular part of the matrix, and the
  lower-triangular part of array A is not referenced.

  When uplo specifies the lower  portion of the matrix,  the leading n by n
  part of the array contains the lower-triangular part of the matrix, and the
  upper-triangular part of array A is not referenced.

  If diag is equal to 'U' or 'u', the diagonal elements of A are also not
  referenced, but are assumed to be unity.
  On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      MAX(1,n).
                      On exit, lda is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is overwritten with the transformed vector
                      x.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

Description  
  The _TRMV subprograms compute a matrix-vector product for a triangular
  matrix or its transpose: x  =  Ax or  x  =  transp(A)*x .  In addition to
  these operations, the CTRMV and ZTRMV subprograms compute a matrix-vector
  product for conjugate transpose: x  =  conjug_transp(A)*x .

  x is a vector with n elements, and A is an n by n, unit or non-unit, upper-
  or lower-triangular matrix.

Example  

  REAL*4 A(50,20), X(20)
  INCX = 1
  N = 20
  LDA = 50
  CALL STRMV('U','N','N',N,A,LDA,X,INCX)

  This FORTRAN code computes the product x  =  Ax where A is an upper-
  triangular matrix, of order 20, with a non-unit diagonal.


**************************************************************************************


strsv, dtrsv, ctrsv, ztrsv 
Solver of a system of linear equations with a triangular matrix

FORMAT  
  {S,D,C,Z}TRSV (uplo, trans, diag, n, a, lda, x, incx)

Arguments  
  uplo                character*1
                      On entry, specifies whether the matrix A is an upper-
                      or lower-triangular matrix:

                      If uplo = 'U' or 'u', A is an upper-triangular matrix.

                      If uplo = 'L' or
                      On exit, uplo is unchanged.

  trans               character*1
                      On entry, specifies the system to be solved:

                      If trans = 'N' or 'n', the system is Ax = b.

                      If trans = 'T' or 't', the system is transp(A)*x = b.

                      If trans = 'C' or 'c', the system is conjug_transp(A)*x
                      = b.
                      On exit, trans is unchanged.

  diag                character*1
                      On entry, specifies whether the matrix A is unit-
                      triangular:

                      If diag = 'U' or 'u', A is a unit-triangular matrix.

                      If diag = 'N' or 'n', A is not a unit-triangular
                      matrix.
                      On exit, diag is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n.

  When uplo specifies the upper portion of the matrix, the leading n by n
  part of the array contains the upper-triangular part of the matrix, and the
  lower-triangular part of array A is not referenced.

  When uplo specifies the lower  portion of the matrix,  the leading n by n
  part of the array contains the lower-triangular part of the matrix, and the
  upper-triangular part of array A is not referenced.

  If diag is equal to 'U' or 'u', the diagonal elements of A are also not
  referenced, but are assumed to be unity.
  On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      MAX(1,n).
                      On exit, lda is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector b.
                      On exit, x is overwritten with the solution vector x.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

Description  
  The _TRSV subprograms solve one of the following systems of linear
  equations for x:  Ax = b  or  transp(A)*x = b .  In addition to these
  operations, the CTRSV and ZTRSV subprograms solve the following systems of
  linear equation:  conjug_transp(A)*x = b .

  b and x are vectors with n elements and A is an n by n, unit or non-unit,
  upper- or lower-triangular matrix.

  The _TRSV routines do not perform checks for singularity or near
  singularity of the triangular matrix.  The requirements for such a test
  depend on the application.  If necessary, perform the test in your
  application program before calling this routine.

Example  

  REAL*8 A(100,40), X(40)
  INCX = 1
  N = 40
  LDA = 100
  CALL DTRSV('L','N','U',N,A,LDA,X,INCX)

  This FORTRAN code solves the system Ax=b where A is a lower-triangular
  matrix of order 40, with a unit diagonal.  The right hand side b is
  originally stored in the vector x.


**************************************************************************************


sger, dger, cgerc, zgerc, cgeru, zgeru 
Rank-one update of a general matrix

FORMAT  
  {S,D}GER (m, n, alpha, x, incx, y, incy, a, lda) {C,Z}GER{C,U} (m, n,
  alpha, x, incx, y, incy, a, lda)

Arguments  
  m                   integer*4
                      On entry, the number of rows of the matrix A; m >= 0.
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns of the matrix A; n >=
                      0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(m-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array of length at least
                      (1+(n-1)*|incy|).  The incremented array Y must contain
                      the vector y.
                      On exit, y is unchanged.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n. The leading m by n part of the array contains the
                      elements of the matrix A.
                      On exit, a is overwritten by the updated matrix.

  lda                 integer*4
                      On entry, the first dimension of A; lda >= MAX(1,m).
                      On exit, lda is unchanged.

Description  
  SGER and DGER perform a rank-one update of a real general matrix: A  =
  alpha*x*transp(y) + A

  CGERU and ZGERU perform a rank-one update of an unconjugated complex
  general matrix: A  = alpha*x*transp(y) + A

  CGERC and ZGERC perform a rank-one update of a conjugated complex general
  matrix: A  = alpha*x*conjug_transp(y) + A

  alphais a scalar, x is an m-element vector, y is an n-element vector, and A
  is an m by n matrix.

Example  

  REAL*4 A(10,10), X(10), Y(5), alpha
  INCX = 1
  INCY = 1
  LDA = 10
  M = 3
  N = 4
  alpha = 2.3
  CALL SGER(M,N,alpha,X,INCX,Y,INCY,A,LDA)

  This FORTRAN code computes the rank-1 update A  =  alpha*x*transp(y)
   + A.  Only the upper left submatrix of A, of dimension (3,4) and starting
  at location A(1,1), is updated.

  COMPLEX*8 A(10,10), X(10), Y(5), alpha
  INCX = 1
  INCY = 1
  LDA = 10
  M = 3
  N = 4
  alpha = (2.3, 1.2)
  CALL CGERC(M,N,alpha,X,INCX,Y,INCY,A,LDA)

  This FORTRAN code computes the rank-1 update A  =  alpha*x*conjug_transp(y)
  + A. Only the upper left submatrix of A, of dimension (3,4) and starting at
  location A(1,1), is updated.


**************************************************************************************


ssyr, dsyr, cher, zher 
Rank-one update of a symmetric or hermitian matrix

FORMAT  
  {S,D}SYR (uplo, n, alpha, x, incx, a, lda) {C,Z}HER (uplo, n, alpha, x,
  incx, a, lda)

Arguments  

  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the array A is referenced:

                      If uplo = 'U' or 'u', the upper-triangular part of A is
                      referenced.

                      If uplo = 'L' or 'l', the lower-triangular part of A is
                      referenced.
                      On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A and the number of
                      elements in vector x; n >= 0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n.

  When uplo specifies the upper portion of the matrix, the leading n by n
  part of the array contains the upper-triangular part of the matrix, and the
  lower-triangular part of array A is not referenced.

  When uplo specifies the lower  portion of the matrix,  the leading n by n
  part of the array contains the lower-triangular part of the matrix, and the
  upper-triangular part of array A is not referenced.

  For CHER and ZHER routines, the imaginary parts of the diagonal elements
  are not accessed, need not be set, and are assumed to be zero.

  On exit, a is overwritten; the specified part of the array A is overwritten
  by the part of the updated matrix.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      MAX(1,n).
                      On exit, lda is unchanged.
Description  SSYR and DSYR perform the rank-one update of a real symmetric matrix: A  =
  alpha*x*transp(x) + A

  CHER and ZHER perform the rank-one update of a complex Hermitian matrix: A
  =  alpha*x*conjug_transp(x) + A

  alpha is a scalar, x is vector with n elements, and A is an n by n matrix
  in packed form. In the case of SSYR and DSYR, matrix A is a symmetric
  matrix and in the case of CHER and ZHER, matrix A is a Hermitian matrix.

Example  

  REAL*4 A(50,20), X(20), alpha
  INCX = 1
  LDA = 50
  N = 20
  alpha = 2.0
  CALL SSYR('L',N,alpha,X,INCX,A,LDA)

  This FORTRAN code computes the rank-1 update of the matrix A, given by A  =
  alpha*x*transp(x)
   + A.  A is a real symmetric matrix with its lower-triangular part stored.

  COMPLEX*16 A(50,20), X(20), alpha
  INCX = 1
  LDA = 50
  N = 20
  alpha = (2.0D0, 1.0D0)
  CALL ZHER('L',N,alpha,X,INCX,A,LDA)

  This FORTRAN code computes the rank-1 update of the matrix A, given by A  =
  alpha*x*conjug_transp(x) + A.  A is a complex Hermitian matrix with its
  lower-triangular part stored.


**************************************************************************************
 

sspr, dspr, chpr, zhpr 
Rank-one update of a symmetric or hermitian matrix stored in packed form

FORMAT  
  {S,D}SPR (uplo, n, alpha, x, incx, ap) {C,Z}HPR (uplo, n, alpha, x, incx,
  ap)

Arguments  
  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the matrix A is supplied in the
                      packed array AP:

                      If uplo = 'U' or 'u', the upper-triangular part of A is
                      supplied.

                      If uplo = 'L' or 'l', the lower-triangular part of A is
                      supplied.
                      On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  ap                  real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array AP of length at least
                      n(n + 1)/2

  If uplo specifies the upper triangular part of the matrix A, the array
  contains those elements of the matrix, packed sequentially, column by
  column, so that AP(1) contains a(11), AP(2) and AP(3) contain a(12) and
  a(22) respectively, and so on.

  If uplo specifies the lower triangular part to the matrix A, the array
  contains those elements of the matrix, also packed sequentially, so that
  AP(1) contains a(11), AP(2) and AP(3) contain a(21) and a(31) respectively,
  and so on.

  For CHPR and ZHPR routines, the imaginary parts of the diagonal elements
  are not accessed, need not be set, and are assumed to be zero.

  On exit, ap is overwritten by the specified part of the updated matrix.

Description  
  SSPR and DSPR perform the rank-one update of a real symmetric matrix stored
  in packed form: A  =  alpha*x*transp(x) + A

  CHPR and ZHPR perform the rank-one update of a complex Hermitian matrix
  stored in packed form: A  =  alpha*x*conjug_transp(x) + A

  alpha is a scalar, x is vector with n elements, and A is an n by n matrix
  in packed form. In the case of SSPR and DSPR, matrix A is a symmetric
  matrix and in the case of CHPR and ZHPR, matrix A is a Hermitian matrix.

Example

  REAL*8 AP(500), X(30), Y(30), alpha
  INCX = 1
  alpha = 1.0D0
  N = 30
  CALL DSPR('U',N,alpha,X,INCX,AP)

  This FORTRAN code computes the rank-1 update A  =  x*transp(x)
   + A where A is a real symmetric matrix, of order 30, with its upper-
  triangular part stored in packed form in AP.


**************************************************************************************
 

ssyr2, dsyr2, cher2, zher2 
Rank-two update of a symmetric or hermitian matrix

FORMAT  
  {S,D}SYR2 (uplo, n, alpha, x, incx, y, incy, a, lda) {C,Z}HER2 (uplo, n,
  alpha, x, incx, y, incy, a, lda)

Arguments  
  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the array A is referenced:

                      If uplo = 'U' or 'u', the upper-triangular part of A is
                      referenced.

                      If uplo = 'L' or 'l', the lower-triangular part of A is
                      referenced.
                      On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).  The incremented array Y must contain
                      the vector y.
                      On exit, y is unchanged.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with dimensions lda
                      by n.

  When uplo specifies the upper portion of the matrix, the leading n by n
  part of the array contains the upper-triangular part of the matrix, and the
  lower-triangular part of array A is not referenced.

  When uplo specifies the lower  portion of the matrix,  the leading n by n
  part of the array contains the lower-triangular part of the matrix, and the
  upper-triangular part of array A is not referenced.

  For complex routines, the imaginary parts of the diagonal elements need not
  be set.  They are assumed to be 0, and on exit they are set to 0.

  On exit, a is overwritten; the specified part of the array A is overwritten
  by the specified part of the updated matrix.

  lda                 integer*4
                      On entry, the first dimension of array A; lda >=
                      MAX(1,n).
                      On exit, lda is unchanged.

Description
  SSYR2 and DSYR2 perform the rank-two update of a real symmetric matrix: A
  =  alpha*x*transp(y)
   + alpha*y*transp(x) + A

  CHER2 and ZHER2 perform the rank-two update of a complex Hermitian matrix:
  A  =  alpha*x*conjug_transp(y) + conjugate(alpha)*y*conjug_transp(x) + A

  alpha is a scalar, x and y are vectors with n elements, and A is an n by n
  matrix. In the case of SSYR2 and DSYR2, matrix A is a symmetric matrix and
  in the case of CHER2 and ZHER2, matrix A is a Hermitian matrix.

Example  

  REAL*8 A(50,20), X(20), Y(20), alpha
  INCX = 1
  LDA = 50
  N = 20
  INCY = 1
  alpha = 1.0D0
  CALL DSYR2('U',N,alpha,X,INCX,Y,INCY,A,LDA)

  This FORTRAN code computes the rank-2 update of a real symmetric matrix A,
  given by A  =  x*transp(y)
   + y*transp(x) + A.  Only the upper-triangular part of A is stored.


**************************************************************************************


sspr2, dspr2, chpr2, zhpr2, 
Rank-two update of a symmetric or hermitian matrix stored in packed form

FORMAT  
  {S,D}SPR2 (uplo, n, alpha, x, incx, y, incy, ap) {C,Z}HPR2 (uplo, n, alpha,
  x, incx, y, incy, ap)

Arguments  
  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the matrix A is supplied in the
                      packed array AP:

                      If uplo = 'U' or 'u', the upper-triangular part of
                      matrix A is supplied.

                      If uplo = 'L' or 'l', the lower-triangular part of
                      matrix A is supplied.
                      On exit, uplo is unchanged.

  n                   integer*4
                      On entry, the order of the matrix A; n >= 0.
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, the scalar alpha*.
                      On exit, alpha is unchanged.

  x                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|).  Array X contains the vector x.
                      On exit, x is unchanged.

  incx                integer*4
                      On entry, the increment for the elements of X; incx
                      must not equal zero.
                      On exit, incx is unchanged.

  y                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).  The incremented array Y must contain
                      the vector y.
                      On exit, y is unchanged.

  incy                integer*4
                      On entry, the increment for the elements of Y; incy
                      must not equal zero.
                      On exit, incy is unchanged.

  ap                  real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array AP of length at least
                      n(n + 1)/2.

  If uplo specifies the upper triangular part of the matrix A, the array
  contains those elements of the matrix, packed sequentially, column by
  column, so that AP(1) contains a(11), AP(2) and AP(3) contain a(12) and
  a(22) respectively, and so on.

  If uplo specifies the lower triangular part to the matrix A, the array
  contains those elements of the matrix, also packed sequentially, so that
  AP(1) contains a(11), AP(2) and AP(3) contain a(21) and a(31) respectively,
  and so on.

  For CHPR2 and ZHPR2 routines, the imaginary parts of the diagonal elements
  are not accessed, need not be set, and are assumed to be zero.

  On exit, ap is overwritten by the specified part of the updated matrix.

Description  
  SSPR2 and DSPR2 perform the rank-two update of a real symmetric matrix
  stored in packed form: A  =  alpha*x*transp(y)
   + alpha*y*transp(x) + A

  CHPR2 and ZHPR2 perform the rank-two update of a complex Hermitian matrix
  stored in packed form: A  =  alpha*x*conjug_transp(y) +
  conjugate(alpha)*y*conjug_transp(x) + A

  alpha is a scalar, x is vector with n elements, and A is an n by n matrix
  in packed form. In the case of SSPR2 and DSPR2, matrix A is a symmetric
  matrix and in the case of CHPR2 and ZHPR2, matrix A is a Hermitian matrix.

Example  

  REAL*4 AP(250), X(20), Y(20), alpha
  INCX = 1
  INCY = 1
  alpha = 2.0
  N = 20
  CALL SSPR2('L',N,alpha,X,INCX,Y,INCY,AP)

  This FORTRAN code computes the rank-2 update of a real symmetric matrix A,
  given by A  =  alpha*x*transp(y)
   + alpha*y*transp(x) + A.  A is a real symmetric matrix, of order 20, with
  its lower-triangular part stored in packed form in AP.


**************************************************************************************

*/
 
#endif
