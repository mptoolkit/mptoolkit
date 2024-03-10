// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/blas3f.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
  blas3f.h

  C++ interface to BLAS level 3
*/

#if !defined(MPTOOLKIT_COMMON_BLAS3F_H)
#define MPTOOLKIT_COMMON_BLAS3F_H

#include "fortran.h"
#include "restrict.h"
#include "trace.h"
#include <complex>

#if defined(BLAS3_TRACE_DETAILED)
#define TRACE_BLAS3(Msg) TRACE(Msg)
#else
#define TRACE_BLAS3(Msg) DUMMY_TRACE(Msg)
#endif


namespace BLAS
{

using namespace Fortran;

// real

void dgemm(char transa, char transb, integer m, integer n, integer k,
           double alpha, double const* a, integer lda, double const* b,
           integer ldb, double beta, double* restrict c, integer ldc);

void dtrsm(char *side, char *uplo, char *transa, char *diag,
           integer *m, integer *n, double *alpha, const double *A, integer *lda,
           const double *B, integer *ldb);

void dtrmm(char *side, char *uplo, char *transa, char *diag,
           integer *m, integer *n, double *alpha, const double *A, integer *lda,
           const double *B, integer *ldb);

void dsymm(char *side, char *uplo, integer *m, integer *n,
           double *alpha, const double *A, integer *lda, const double *B,
           integer *ldb, double *beta, double *C, integer *ldc);

void dsyrk(char *uplo, char *transa, integer *n, integer *k,
           double *alpha, double *A, integer *lda, double *beta, double *C,
           integer *ldc);

void dsyr2k(char *uplo, char *transa, integer *n, integer *k,
            double *alpha, double *A, integer *lda, double *B, integer *ldb,
            double *beta, double *C, integer *ldc);


// complex

void zgemm(char transa, char transb, integer m, integer n, integer k,
           std::complex<double> alpha, std::complex<double> const* a, integer lda,
           std::complex<double> const* b, integer ldb, std::complex<double> beta,
           std::complex<double>* restrict c, integer ldc);

// the raw functions are in their own namespace, the wrapper functions call these.
namespace raw
{
using Fortran::complex;
extern "C"
{
   void F77NAME(dgemm)(char const* transa, char const* transb, integer const* m, integer const* n,
                       integer const* k,
                       double const* alpha, double const* a, integer const* lda, double const* b,
                       integer const* ldb, double const* beta, double* restrict c, integer const* ldc);

   void F77NAME(dtrsm)(char *side, char *uplo, char *transa, char *diag,
                       integer *m, integer *n, double *alpha, const double *A, integer *lda,
                       const double *B, integer *ldb);

   void F77NAME(dtrmm)(char *side, char *uplo, char *transa, char *diag,
                       integer *m, integer *n, double *alpha, const double *A, integer *lda,
                       const double *B, integer *ldb);

   void F77NAME(dsymm)(char *side, char *uplo, integer *m, integer *n,
                       double *alpha, const double *A, integer *lda, const double *B,
                       integer *ldb, double *beta, double *C, integer *ldc);

   void F77NAME(dsyrk)(char *uplo, char *transa, integer *n, integer *k,
                       double *alpha, double *A, integer *lda, double *beta, double *C,
                       integer *ldc);

   void F77NAME(dsyr2k)(char *uplo, char *transa, integer *n, integer *k,
                        double *alpha, double *A, integer *lda, double *B, integer *ldb,
                        double *beta, double *C, integer *ldc);

   // complex

   void F77NAME(zgemm)(char const* transa, char const* transb,
                       integer const* m, integer const* n, integer const* k,
                       std::complex<double> const* alpha, std::complex<double> const* a, integer const* lda,
                       std::complex<double> const* b, integer const* ldb, std::complex<double> const* beta,
                       std::complex<double>* restrict c, integer const* ldc);


} // extern "C"

} // namespace raw

#if defined(DEBUG_DGEMM)
namespace DGS
{
   extern char transa;
   extern char transb;
   extern integer m;
   extern integer n;
   extern integer k;
   extern double alpha;
   extern double const* a;
   extern integer lda;
   extern double const* b;
   extern integer ldb;
   extern double beta;
   extern double const* c;
   extern integer ldc;

   void DebugPrintDgemm();
} // namespace DGS
#endif

// inlines

// real

inline void dgemm(char transa, char transb, integer m, integer n, integer k,
                  double alpha, double const* a, integer lda, double const* b,
                  integer ldb, double beta, double* restrict c, integer ldc)
{
   TRACE_BLAS3("BLAS3: dgemm")(transa)(transb)(m)(n)(k)(alpha)(a)(lda)(b)(ldb)(beta)(c)(ldc);
#if defined(DEBUG_DGEMM)
   DGS::transa = transa;
   DGS::transb = transb;
   DGS::m = m;
   DGS::n = n;
   DGS::k = k;
   DGS::alpha = alpha;
   DGS::a = a;
   DGS::lda = lda;
   DGS::b = b;
   DGS::ldb = ldb;
   DGS::beta = beta;
   DGS::c = c;
   DGS::ldc = ldc;
#endif

    raw::F77NAME(dgemm)(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
}

inline void dtrsm(char *side, char *uplo, char *transa, char *diag,
                  integer *m, integer *n, double *alpha, const double *A, integer *lda,
                  const double *B, integer *ldb)
{
    raw::F77NAME(dtrsm)(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);
}

inline void dtrmm(char *side, char *uplo, char *transa, char *diag,
                  integer *m, integer *n, double *alpha, const double *A, integer *lda,
                  const double *B, integer *ldb)
{
   raw::F77NAME(dtrmm)(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);
}

inline void dsymm(char *side, char *uplo, integer *m, integer *n,
                  double *alpha, const double *A, integer *lda, const double *B,
                  integer *ldb, double *beta, double *C, integer *ldc)
{
   raw::F77NAME(dsymm)(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

inline void dsyrk(char *uplo, char *transa, integer *n, integer *k,
                  double *alpha, double *A, integer *lda, double *beta, double *C,
                  integer *ldc)
{
   raw::F77NAME(dsyrk)(uplo, transa, n, k, alpha, A, lda, beta, C, ldc);
}

inline void dsyr2k(char *uplo, char *transa, integer *n, integer *k,
                   double *alpha, double *A, integer *lda, double *B, integer *ldb,
                   double *beta, double *C, integer *ldc)
{
   raw::F77NAME(dsyr2k)(uplo, transa, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

// complex

inline
void zgemm(char transa, char transb, integer m, integer n, integer k,
           std::complex<double> alpha, std::complex<double> const* a, integer lda,
           std::complex<double> const* b, integer ldb, std::complex<double> beta,
           std::complex<double>* restrict c, integer ldc)
{
   //   using Fortran::complex;
   TRACE_BLAS3("BLAS3: zgemm")(transa)(transb)(m)(n)(k)(alpha)(a)(lda)(b)(ldb)(beta)(c)(ldc);
   raw::F77NAME(zgemm)(&transa, &transb, &m, &n, &k,
                       &alpha, a, &lda, b, &ldb, &beta, c, &ldc);
   //                  reinterpret_cast<complex const*>(&alpha),
   //                  reinterpret_cast<complex const*>(a), &lda,
   //                  reinterpret_cast<complex const*>(b), &ldb,
   //                  reinterpret_cast<complex const*>(&beta),
   //                  reinterpret_cast<complex*>(c), &ldc);
}

} // namespace BLAS


/*
BLAS documentation

**************************************************************************************


sgemm, dgemm, cgemm, zgemm
Matrix-matrix product and addition

FORMAT
  {S,D,C,Z}GEMM (transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

Arguments

  transa              character*1
                      On entry, specifies the form of (op)A used in the
                      matrix multiplication:

                      If transa = 'N' or 'n', (op)A = A

                      If transa = 'T' or 't', (op)A = transp(A)

                      If transa = 'R' or 'r', (op)A = conjugate(A)

                      If transa = 'C' or 'c', (op)A = conjug_transp(A)
                      On exit, transa is unchanged.

  transb              character*1
                      On entry, specifies the form of (op)B used in the
                      matrix multiplication:

                      If transb = 'N' or 'n', (op)B = B

                      If transb = 'T' or 't', (op)B = transp(B)

                      If transb = 'R' or 'r', (op)B = conjugate(B)

                      If transb = 'C' or 'c', (op)B = conjug_transp(B)

  m                   integer*4
                      On entry, the number of rows of the matrix (op)A and of
                      the matrix C; m >= 0
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns of the matrix (op)B and
                      of the matrix C; n >= 0
                      On exit, n is unchanged.

  k                   integer*4
                      On entry, the number of columns of the matrix (op)A and
                      the number of rows of the matrix (op)B; k >= 0
                      On exit, k is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar alpha.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array A with dimensions lda
                      by ka.
                      For (op)A = A  or  conjugate(A), ka >= k and the lead-
                      ing m by k portion of the array A contains the matrix
                      A.
                      For (op)A = transp(A)   or conjug_transp(A), ka >= m
                      and the leading k by m part of the array A contains the
                      matrix A.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A.
                      For (op)A = A  or conjugate(A), lda >= MAX(1,m).
                      For (op)A = transp(A)  or conjug_transp(A), lda >=
                      MAX(1,k).
                      On exit, lda is unchanged.

  b                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array B with dimensions ldb
                      by kb.
                      For (op)B = B or conjugate(B), kb >= n and the leading
                      k by n portion of the array contains the matrix B.
                      For (op)B = transp(B) or conjug_transp(B), kb >= k and
                      the leading n by k part of the array contains the
                      matrix B.
                      On exit, b is unchanged.

  ldb                 integer*4
                      On entry, the first dimension of array B.
                      For (op)B = B or = MAX(1,k).
                      For (op)B = transp(B) or conjug_transp(B), ldb >=
                      MAX(1,n).
                      On exit, ldb is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar beta.
                      On exit, beta is unchanged.

  c                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with the dimension
                      ldc by at least n.
                      On exit,  the leading  m by n part of array C is
                      overwritten by the matrix alpha*(op)A*(op)B + beta*C.

  ldc                 integer*4
                      On entry, the first dimension  of array C; ldc >=
                      MAX(1,m)
                      On exit, ldc is unchanged.

Description
  The _GEMM routines perform the following operations: C  = alpha(op)A(op)B +
  beta*C
  where (op)(X) = X, transp(X), conjugate(X),  or conjug_transp(X), alpha and
  beta are scalars, and A, B, and C are matrices. (op)A is an m by k matrix,
  (op)B is a k by n matrix, and C is an m by n matrix.

Example

  REAL*4 A(20,40), B(20,30), C(40,30), alpha, beta
  M = 10
  N = 20
  K = 15
  LDA = 20
  LDB = 20
  LDC = 40
  alpha = 2.0
  beta = 2.0
  CALL SGEMM ('T','N',M,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)

  This FORTRAN code computes the product C  =  alpha * transp(A)*B + beta*C
  where A is a real general matrix.  A is a 15 by 10 real general matrix
  embedded in array A.  B is a 15 by 20 real general matrix embedded in array
  B.  C is a  10 by 20 real general matrix embedded in array C.


**************************************************************************************


strsm, dtrsm, ctrsm, ztrsm
Solve a triangular system of equations with a triangular coefficient matrix

FORMAT
  {S,D,C,Z}TRSM ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )

Arguments

  side                character*1
                      On entry, specifies whether (op)A is on the left side
                      or the right side of X in the system of equations:

                      If side = 'L' or 'l', the system is (op)A X   =   alpha
                      * B.

                      If side = 'R' or 'r', the system is X (op)A   =  alpha
                      * B.
                      On exit, side is unchanged.

  uplo                character*1
                      On entry, specifies whether the matrix A is an upper-
                      or lower-triangular matrix:

                      If uplo = 'U' or 'u', the matrix A is an upper-
                      triangular matrix.

                      If uplo = 'L' or 'l', the matrix A is a lower-
                      triangular matrix.
                      On exit, uplo is unchanged.

  transa              character*1
                      On entry, specifies the form of (op)A used in the sys-
                      tem of equations:

                      If transa = 'N' or 'n', (op)A = A.

                      If transa = 'T' or 't', (op)A = transp(A).

                      If transa = 'C' or 'c', (op)A = conjug_transp(A).
                      On exit, transa is unchanged.

  diag                character*1
                      On entry, specifies whether the matrix A is unit-
                      triangular:

                      If diag = 'U' or 'u', A is a unit-triangular matrix.

                      If diag = 'N' or 'n', A is not a unit-triangular
                      matrix.
                      On exit, diag is unchanged.

  m                   integer*4
                      On entry, the number of rows m of the matrix B; m >= 0
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns n of the matrix B; n >=
                      0
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar alpha.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array A with dimensions lda
                      by k.
                      If the multiplication is on the left side, k >= m and
                      the leading m by m part of the array contains the
                      matrix A.
                      If the multiplication is on the right side, k >= n and
                      the leading n by n part of the array A must contain the
                      matrix A.
                      In either case, when the leading part of the array is
                      specified as the  upper part, the upper triangular part
                      of array A contains the upper-triangular part of the
                      matrix A, and the lower-triangular part of matrix A is
                      not referenced.  When the lower part is specified, the
                      lower triangular part of the array A contains the lower
                      triangular part of the matrix  A, and the upper-
                      triangular part of A is not  referenced.

                      If matrix A is unit-triangular, its diagonal elements
                      are assumed to be unity and are not referenced.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of A.  When multiplica-
                      tion is on the left, lda >= MAX(1,m). When multiplica-
                      tion is on the right, lda >= MAX(1,n).
                      On exit, lda is unchanged.

  b                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array B of dimensions ldb
                      by at least n.  The leading m by n part of the array B
                      must contain the right-hand-side matrix B.
                      On exit, b is overwritten by the m by n solution matrix
                      X.

  ldb                 integer*4
                      On entry, the first dimension of B; ldb >= MAX(1,m)
                      On exit, ldb is unchanged.

Description
  The _TRSM routines solve a triangular system of equations where the coeffi-
  cient matrix A is a triangular matrix: (op)AX  =  alpha * B X(op)A   =
  alpha * B
  (op)A = A, transp(A),  or  conjug_transp(A) , alpha is a scalar, X and B
  are m by n matrices, and A is a unit or non-unit, upper- or lower-
  triangular matrix.

Example

  REAL*8 A(100,40), B(40,20), alpha
  M = 16
  N = 18
  LDA = 100
  LDB = 40
  alpha = 2.0D0
  CALL DTRSM ('L','U','N','U',M,N,alpha,A,LDA,B,LDB)

  This FORTRAN code solves the system AX=alpha * B where A is an upper-
  triangular real matrix with a unit diagonal.  X and B are 16 by 18
  matrices.  The leading 16 by 16 upper-triangular part of the array A must
  contain the upper-triangular matrix A.  The leading 16 by 18 part of the
  array B must contain the matrix B.  The lower-triangular part of A and the
  diagonal are not referenced.  The leading 16 by 18 part of B is overwritten
  by the solution matrix X.



**************************************************************************************


strmm, dtrmm, ctrmm, ztrmm
Matrix-matrix product for triangular matrix

FORMAT
  {S,D,C,Z}TRMM ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )

Arguments

  side                character*1
                      On entry, specifies whether (op)A multiplies B on the
                      left or right in the operation:

                      If side = 'L' or 'l', the operation is B  =  alpha *
                      (op)A*B.

                      If side = 'R' or 'r', the operation is B  =  alpha * B
                      * (op)A .
                      On exit, side is unchanged.

  uplo                character*1
                      On entry, specifies whether the matrix A is an upper-
                      or lower-triangular matrix:

                      If uplo = 'U' or 'u', the matrix A is an upper-
                      triangular matrix.

                      If uplo = 'L' or 'l', the matrix A is a lower-
                      triangular matrix.
                      On exit, uplo is unchanged.

  transa              character*1
                      On entry, specifies the form of (op)A used in the
                      matrix multiplication:

                      If transa = 'N' or 'n', (op)A = A.

                      If transa = 'T' or 't', (op)A = transp(A).

                      If transa = 'C' or 'c', (op)A = conjug_transp(A).
                      On exit, transa is unchanged.

  diag                character*1
                      On entry, specifies whether the matrix A is unit-
                      triangular:

                      If diag = 'U' or 'u', A is a unit-triangular matrix.

                      If diag = 'N' or 'n', A is not a unit-triangular
                      matrix.
                      On exit, diag is unchanged.

  m                   integer*4
                      On entry, the number of rows of the matrix B; m >= 0
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns of the matrix B; n >= 0
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar alpha.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array A with dimensions lda
                      by k.
                      If the multiplication is on the left side, k >= m and
                      the leading m by m part of the array contains the
                      matrix A.
                      If the multiplication is on the right side, k >= n and
                      the leading n by n part of the array A must contain the
                      matrix A.
                      In either case, when the leading part of the array is
                      specified as the  upper part, the upper triangular part
                      of array A contains the upper-triangular part of the
                      matrix A, and the lower-triangular part of matrix A is
                      not referenced.  When the lower part is specified, the
                      lower triangular part of the array A contains the lower
                      triangular part of the matrix  A, and the upper-
                      triangular part of A is not  referenced.

                      If matrix A is unit-triangular, its diagonal elements
                      are assumed to be unity and are not referenced.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of A.  When multiplica-
                      tion is on the left, lda >= MAX(1,m). When multiplica-
                      tion is on the right, lda >= MAX(1,n).
                      On exit, lda is unchanged.

  b                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array B of dimensions ldb
                      by at least n.  The leading m by n part of the array B
                      must contain the matrix B.
                      On exit, b is overwritten by the m by n updated matrix.

  ldb                 integer*4
                      On entry, the first dimension of B; ldb >= MAX(1,m)
                      On exit, ldb is unchanged.

Description
  STRMM and DTRMM compute a matrix-matrix product for a real triangular
  matrix or its transpose.  CTRMM and ZTRMM compute a matrix-matrix product
  for a complex triangular matrix, its transpose, or its conjugate transpose.
    B = alpha(op)A*B
    B = alpha * B((op)A)

  where (op)A = A, transp(A),  or conjug_transp(A)

  alpha is a scalar, B is an m by n matrix, and A is a unit or non-unit,
  upper- or lower-triangular matrix.

Example

  REAL*8 A(25,40), B(30,35), alpha
  M = 15
  N = 18
  LDA = 25
  LDB = 30
  alpha = -1.0D0
  CALL DTRMM ('R','L','T','U',M,N,alpha,A,LDA,B,LDB)

  This FORTRAN code solves the system B  =  alpha * B*transp(A) where A is a
  lower-triangular real matrix with a unit diagonal.  A is an 18 by 18 real
  triangular matrix embedded in array A, and B is a 15 by 18 real rectangular
  matrix embedded in array B.  The leading 18 by 18 lower-triangular part of
  the array A must contain the lower-triangular matrix A.  The upper-
  triangular part of A and the diagonal are not referenced.

  COMPLEX*16 A(25,40), B(30,35), alpha
  M = 15
  N = 18
  LDA = 25
  LDB = 30
  alpha = (-1.0D0, 2.0D0)
  CALL ZTRMM ('R','L','T','U',M,N,alpha,A,LDA,B,LDB)

  This FORTRAN code solves the system B  =  alpha * B*transp(A) where A is a
  lower-triangular complex matrix with a unit diagonal.  A is an 18 by 18
  complex triangular matrix embedded in array A, and B is a 15 by 18 complex
  rectangular matrix embedded in array B.  The leading 18 by 18 lower-
  triangular part of the array A must contain the lower-triangular matrix A.
  The upper-triangular part of A and the diagonal are not referenced.


**************************************************************************************


ssymm, dsymm,   csymm, zsymm, chemm, zhemm - Matrix-matrix product and addition for a symmetric or hermitian matrix

FORMAT
  {S,D,C,Z}SYMM ( side, uplo,   m, n, alpha, a, lda, b, ldb, beta, c, ldc )
  {C,Z}HEMM ( side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc )

Arguments

  side                character*1
                      On entry, specifies whether the symmetric matrix A mul-
                      tiplies B on the left side or the right side:

                      If side = 'L' or 'l', the operation is C  =  alpha *
                      A*B + beta*C.

                      If side = 'R' or 'r', the operation is C  =  alpha *
                      B*A + beta*C.
                      On exit, side is unchanged.

  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the symmetric matrix A is refer-
                      enced:

                      If uplo = 'U' or 'u', the upper-triangular part of A is
                      referenced.

                      If uplo = 'L' or 'l', the lower-triangular part of A is
                      referenced.
                      On exit, uplo is unchanged.

  m                   integer*4
                      On entry, the number of rows of the matrix C; m >= 0
                      On exit, m is unchanged.

  n                   integer*4
                      On entry, the number of columns of the matrix C; n >= 0
                      On exit, n is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar alpha.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array A with dimensions lda
                      by ka.
                      If the multiplication is on the left side, ka >= m and
                      the leading m by m part of the array contains the
                      matrix A.
                      If the multiplication is on the right side, ka >= n and
                      the leading n by n part of the array A must contain the
                      matrix A.
                      In either case, when the leading part of the array is
                      specified as the  upper part, the upper triangular part
                      of array A contains the upper-triangular part of the
                      matrix A, and the lower-triangular part of matrix A is
                      not referenced.  When the lower part is specified, the
                      lower triangular part of the array A contains the lower
                      triangular part of the matrix  A, and the upper-
                      triangular part of A is not  referenced.

                      In complex Hermitian matrices, the imaginary parts of
                      the diagonal elements need not be set.  They are
                      assumed to be zero.

                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A.  When multi-
                      plication is on the left, lda >= MAX(1,m). When multi-
                      plication is on the right, lda >= MAX(1,n).
                      On exit, lda is unchanged.

  b                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array B of dimensions ldb
                      by at least n.  The leading m by n part of the array B
                      must contain the matrix B.
                      On exit, b is unchanged.

  ldb                 integer*4
                      On entry, the first dimension of B; ldb >= MAX(1,m)
                      On exit, ldb is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar beta.
                      On exit, beta is unchanged.

  c                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array with the dimension
                      ldc by at least n.
                      On exit, c is overwritten; the array C is overwritten
                      by the m by n updated matrix.

  ldc                 integer*4
                      On entry, the first dimension  of array C; ldc >=
                      MAX(1,n)
                      On exit, ldc is unchanged.

Description
  These routines compute a matrix-matrix product and addition for a real or
  complex symmetric matrix or a complex Hermitian matrix:
  C  = alpha * A*B + beta*C
  alpha and beta are scalars, A is the symmetric or Hermitian matrix, and B
  and C are m by n matrices.

Example

  REAL*4 A(20,20), B(30,40), C(30,50), alpha, beta
  M = 10
  N = 20
  LDA = 20
  LDB = 30
  LDC = 30
  alpha = 2.0
  beta = 3.0
  CALL SSYMM ('L','U',M,N,alpha,A,LDA,B,LDB,beta,C,LDC)

  This FORTRAN code computes the product of a symmetric matrix and a rec-
  tangular matrix. The operation is C  =  alpha * A*B + beta*C where A is a
  10 by 10 real symmetric matrix embedded in array A, B is a 10 by 20 real
  matrix embedded in array B, and C is a  10 by 20 real matrix embedded in
  array C.  The leading 10 by 10 upper-triangular part of the array A con-
  tains the upper-triangular part of the matrix A.  The lower-triangular part
  of A is not referenced.

  COMPLEX*16 A(30,40), B(15,20), C(19,13), alpha, beta
  M = 12
  N = 7
  LDA = 30
  LDB = 15
  LDC = 19
  alpha = (2.0D0, 0.0D0)
  beta = (0.0D0, -2.0D0)
  CALL ZHEMM ('R','L',M,N,alpha,A,LDA,B,LDB,beta,C,LDC)

  This FORTRAN code computes the product of a Hermitian matrix and a rec-
  tangular matrix. The operation is C  =  alpha * B*A + beta*C where A is a 7
  by 7 complex Hermitian matrix embedded in array A, B is a 12 by 7 complex
  matrix embedded in array B, and C is a  12 by 7 complex matrix embedded in
  array C.  The leading 7 by 7 lower-triangular part of the array A contains
  the lower-triangular part of the matrix A.  The upper-triangular part of A
  is not referenced.


**************************************************************************************


ssyrk, dsyrk, csyrk, zsyrk
Rank-k update of a symmetric matrix

FORMAT
  {S,D,C,Z}SYRK ( uplo, trans, n, k, alpha, a, lda, beta, c, ldc )

Arguments

  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the symmetric matrix C is to be
                      referenced:

                      If uplo = 'U' or 'u', the upper-triangular part of C is
                      to be referenced.

                      If uplo = 'L' or 'l', the lower-triangular part of C is
                      to be referenced.
                      On exit, uplo is unchanged.

  trans               character*1
                      On entry, specifies the operation to be performed:

                      If trans = 'N' or 'n', C  =  alpha * A*transp(A) +
                      beta*C

                      If trans = 'T' or 't', C  =  alpha * transp(A)A +
                      beta*C
                      On exit, trans is unchanged.

  n                   integer*4
                      On entry, specifies the order of the matrix C; n >= 0
                      On exit, n is unchanged.

  k                   integer*4
                      On entry,  the number of columns of the matrix A when
                      trans = 'N' or the number of rows of the matrix A when
                      trans = 'T' or
                      On exit, k is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar alpha.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array A with dimensions lda
                      by ka.
                      For trans = 'N' or leading n by k portion of the array
                      A contains the matrix A.
                      For trans = 'T' or ka >= n and the leading k by n part
                      of the array A contains the matrix A.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A.
                      For trans = 'N' or lda >= MAX(1,n).
                      For trans = 'T', lda >= MAX(1,k).
                      On exit, lda is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar beta.
                      On exit, beta is unchanged.

  c                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array C of dimensions ldc
                      by at least n.  If uplo specifies the upper part, the
                      leading n by n upper-triangular part of the array C
                      must contain the upper-triangular part of the symmetric
                      matrix C, and the strictly lower-triangular part of C
                      is not referenced.

  If uplo specifies the lower part, the leading n by n lower-triangular part
  of the array C must contain the lower-triangular part of the symmetric
  matrix C, and the strictly upper-triangular part of C is not referenced.
  On exit, c is overwritten; the triangular part of the array C is overwrit-
  ten by the triangular part of the updated matrix.

  ldc                 integer*4
                      On entry, the first dimension  of array C; ldc >=
                      MAX(1,n)
                      On exit, ldc is unchanged.

Description
  The _SYRK routines perform the rank-k update of a symmetric matrix: C  =
  alpha * A*transp(A) + beta*C
  alpha and beta are scalars,  C is an n by n symmetric matrix. In the first
  case, A is an n by k matrix, and in the second case, A is a k by n matrix.

Example

  REAL*4 A(40,20), C(20,20), alpha, beta
  LDA = 40
  LDC = 20
  N = 10
  K = 15
  alpha = 1.0
  beta = 2.0
  CALL SSYRK ('U','N',N,K,alpha,A,LDA,beta,C,LDC)

  This FORTRAN code computes the rank-k update of the real symmetric matrix
  C: C  =  alpha * A*transp(A) + beta*C.  C is a 10 by 10 matrix, and A is a
  10 by 15 matrix.  Only the upper-triangular part of C is referenced.  The
  leading 10 by 15 part of array A contains the matrix A.  The leading 10 by
  10 upper-triangular part of array C contains the upper-triangular matrix C.


**************************************************************************************


ssyr2k, dsyr2k, csyr2k, zsyr2k
Rank-2k update of a symmetric matrix

FORMAT
  {S,D,C,Z}SYR2K ( uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

Arguments

  uplo                character*1
                      On entry, specifies whether the upper- or lower-
                      triangular part of the symmetric matrix C is to be
                      referenced:

                      If uplo = 'U' or 'u', the upper-triangular part of C is
                      to be referenced.

                      If uplo = 'L' or 'l', the lower-triangular part of C is
                      to be referenced.
                      On exit, uplo is unchanged.

  trans               character*1
                      On entry, specifies the operation to be performed:

                      If trans = 'N' or 'n', C  =  alpha * A*transp(B) +
                      alpha * B*transp(A) + beta*C

                      If trans = 'T' or 't', C  =  alpha * transp(A)*B +
                      alpha * transp(B)A + beta*C
                      On exit, trans is unchanged.

  n                   integer*4
                      On entry, the order n of the matrix C; n >= 0
                      On exit, n is unchanged.

  k                   integer*4
                      On entry,  the number of columns of the matrices A and
                      B when trans = 'N' or the number of rows of the matrix
                      A and B when trans = 'T' or k >= 0.
                      On exit, k is unchanged.

  alpha               real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar alpha.
                      On exit, alpha is unchanged.

  a                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array A with dimensions lda
                      by ka.
                      For trans = 'N' or ka >= k and the leading n by k por-
                      tion of the array A contains the matrix A.
                      For trans = 'T' or ka >= n and the leading k by n part
                      of the array A contains the matrix A.
                      On exit, a is unchanged.

  lda                 integer*4
                      On entry, the first dimension of array A.
                      For trans = 'N' or 'n' lda >= MAX(1,n).
                      For trans = 'T' or lda >= MAX(1,k).
                      On exit, lda is unchanged.

  b                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array B with dimensions ldb
                      by kb.
                      For trans = 'N' or kb >= k and the leading n by k por-
                      tion of the array B contains the matrix B.
                      For trans = 'T' or kb >= n and the leading k by n part
                      of the array B contains the matrix B.
                      On exit, b is unchanged.

  ldb                 integer*4
                      On entry, the first dimension of array B.
                      For trans = 'N' or ldb >= MAX(1,n).
                      For trans = 'T' or ldb >= MAX(1,k).
                      On exit, ldb is unchanged.

  beta                real*4 | real*8 | complex*8 | complex*16
                      On entry, specifies the scalar beta.
                      On exit, beta is unchanged.

  c                   real*4 | real*8 | complex*8 | complex*16
                      On entry, a two-dimensional array C of dimensions ldc
                      by at least n.

  If uplo specifies the upper part, the leading n by n upper-triangular part
  of the array C must contain the upper-triangular part of the symmetric
  matrix C, and the strictly lower-triangular part of C is not referenced.

  If uplo specifies the lower part, the leading n by n lower-triangular part
  of the array C must contain the lower-triangular part of the symmetric
  matrix C, and the strictly upper-triangular part of C is not referenced.
  On exit, c is overwritten; the triangular part of the array C is overwrit-
  ten by the triangular part of the updated matrix.

  ldc                 integer*4
                      On entry, the first dimension  of array C; ldc >=
                      MAX(1,n)
                      On exit, ldc is unchanged.

Description
  The _SYR2K routines perform the rank-2k update of a symmetric matrix: C  =
  alpha * A*transp(B) + alpha * B*transp(A)
   + beta*C C  = alpha * transp(A)*B + alpha * transp(B)A
   + beta*C
  alpha and beta are scalars,  C is an n by n symmetric matrix, and A and B
  are n by k matrices in the first case and k by n matrices in the second
  case.

Example

  REAL*4 A(40,10), B(40,10), C(20,20), alpha, beta
  LDA = 40
  LDB = 30
  LDC = 20
  N = 18
  K = 10
  alpha = 1.0
  beta = 2.0
  CALL SSYR2K ('U','N',N,K,alpha,A,LDA,B,LDB,beta,C,LDC)

  This FORTRAN code computes the rank-2k update of the real symmetric matrix
  C: C  =  alpha * A*transp(B) + alpha * B*transp(A) + beta*C.  Only the
  upper-triangular part of C is referenced.  The leading 18 by 10 part of
  array A contains the matrix A.  The leading 18 by 10 part of array B con-
  tains the matrix B.  The leading 18 by 18 upper-triangular part of array C
  contains the upper-triangular matrix C.


**************************************************************************************

*/

#endif
