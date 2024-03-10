// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/eigen.cpp
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
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
//
// implementation of eigen.h LAPACK wrappers
//
// Created 2004-06-07 Ian McCulloch
//

#include "eigen.h"
#include "diagonalmatrix.h"
#include "common/lapackf.h"
#include <algorithm>

namespace LinearAlgebra
{

namespace detail
{
   // put the ARPACK mutex here, avoids having another .cpp file
   std::mutex ArpackMutex;
}

namespace Private
{

void LinearSolveSPD(int Size, int Nrhs, double* A, int ldA, double* B, int ldB)
{
   char uplo = 'L';
   Fortran::integer info = 0;

   // do the actual call
   LAPACK::dposv(uplo, Size, Nrhs, A, ldA, B, ldB, info);

   CHECK(info == 0)("LAPACK::dposv")(info);
}

void LinearSolveHPD(int Size, int Nrhs,
                    std::complex<double>* A, int ldA,
                    std::complex<double>* B, int ldB)
{
   char uplo = 'L';
   Fortran::integer info = 0;

   // do the actual call
   LAPACK::zposv(uplo, Size, Nrhs, A, ldA, B, ldB, info);

   CHECK(info == 0)("LAPACK::zposv")(info);
}

void LinearSolve(int Size, int Nrhs, double* A, int ldA, double* B, int ldB)
{
   Fortran::integer info = 0;

   Fortran::integer* IPiv = new Fortran::integer[Size];
   // do the actual call
   LAPACK::dgesv(Size, Nrhs, A, ldA, IPiv, B, ldB, info);
   delete[] IPiv;

   CHECK(info == 0)("LAPACK::dgesv")(info);
}

void LeastSquares(int M, int N, int Nrhs, double* A, int ldA, double* B, int ldB)
{
   Fortran::integer info = 0;

   Fortran::integer LWork = std::max(1, std::min(M,N) + std::max(std::min(M,N), Nrhs));
   double* Work = new double[LWork];
   // do the actual call
   LAPACK::dgels('N', M, N, Nrhs, A, ldA, B, ldB, Work, LWork, info);
   delete[] Work;

   CHECK(info == 0)("LAPACK::dgels")(info);
}

void EigenvaluesSymmetric(int Size, double* Data, int LeadingDim, double* Eigen)
{
   // FIXME: we can do better than this by calling LAPACK such that we don't
   // bother calculating the eigenvectors
   DiagonalizeSymmetric(Size, Data, LeadingDim, Eigen);
}

void EigenvaluesHermitian(int Size, std::complex<double>* Data,
                          int LeadingDim, double* Eigen)
{
   // FIXME: we can do better than this by calling LAPACK such that we don't
   // bother calculating the eigenvectors
   DiagonalizeHermitian(Size, Data, LeadingDim, Eigen);
}

void EigenvaluesComplex(int Size, std::complex<double>* Data,
                        int LeadingDim, std::complex<double>* Eigen)
{
   // First step: balance the matrix
   Fortran::integer ilo = 0, ihi = 0, info = 0;
   double* Scale = new double[Size];
   LAPACK::zgebal('B', Size, Data, LeadingDim, ilo, ihi, Scale, info);
   CHECK(info == 0)("LAPACK::zgebal")(info);

   // Second step: reduce matrix to upper Hessenberg form
   // firstly determine the optimal size of the workspace
   std::complex<double>* Tau = new std::complex<double>[Size-1];
   std::complex<double> WorkSize;
   Fortran::integer LWork = -1;
   LAPACK::zgehrd(Size, ilo, ihi, Data, LeadingDim, Tau, &WorkSize, -1, info);
   LWork = Fortran::integer(WorkSize.real());
   std::complex<double>* Work = new std::complex<double>[LWork];
   // do the actual call
   LAPACK::zgehrd(Size, ilo, ihi, Data, LeadingDim, Tau, Work, LWork, info);
   CHECK(info == 0)("LAPACK::zgehrd")(info);
   delete[] Work;

   // Third step: compute the eigenvalues of the Hessenberg matrix
   // workspace query
   std::complex<double> z;
   LAPACK::zhseqr('E', 'N', Size, ilo, ihi, Data, LeadingDim, Eigen, &z, 1, &WorkSize, -1, info);
   LWork = Fortran::integer(WorkSize.real());
   Work =  new std::complex<double>[LWork];
   LAPACK::zhseqr('E', 'N', Size, ilo, ihi, Data, LeadingDim, Eigen, &z, 1, Work, LWork, info);
   CHECK(info == 0)("LAPACK::zgesqr")(info);
   delete[] Work;

   delete[] Tau;
   delete[] Scale;
}

// for debugging
Matrix<double> makeMatrix(double const* Data, int Size, int LeadingDim)
{
   Matrix<double> M(Size, Size);
   for (int i = 0; i < Size; ++i)
   {
      for (int j = 0; j < Size; ++j)
      {
         M(i,j) = Data[i*LeadingDim + j];
      }
   }
   return M;
}

void DiagonalizeSymmetric(int Size, double* Data, int LeadingDim, double* Eigen)
{
   char jobz = 'V';
   char uplo = 'L';
   double worksize;
   double* work = &worksize;
   int lwork = -1;                         // for query of lwork
   Fortran::integer info = 0;

   // query for the optimial size of the workspace
   LAPACK::dsyev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, info);

   lwork = int(work[0]);
   work = new double[lwork];

   // do the actual call
   LAPACK::dsyev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, info);

   CHECK(info == 0)("LAPACK::dsyev")(info)(makeMatrix(Data, Size, LeadingDim));

   delete[] work;
}

void DiagonalizeHermitian(int Size, std::complex<double>* Data, int LeadingDim, double* Eigen)
{
   if (Size == 0)
   {
      DEBUG_WARNING("Zero size DiagonalizeHermitian");
      return;
   }

   char jobz = 'V';
   char uplo = 'L';
   std::complex<double>  worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;                         // for query of lwork
   double* rwork;
   int lrwork = std::max(1, 3*Size-2);
   Fortran::integer info = 0;

   rwork = new double[lrwork];

   // query for the optimial size of the workspace
   LAPACK::zheev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zheev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, rwork, info);

   CHECK(info == 0)("LAPACK::zheev")(info);

   delete[] work;
   delete[] rwork;
}

void Diagonalize(int Size, std::complex<double> const* DataX, int LeadingDim, std::complex<double>* Eigen,
                 std::complex<double>* LeftVectors, int LeftLeadingDim,
                 std::complex<double>* RightVectors, int RightLeadingDim)
{
   char jobl = 'V';
   char jobr = 'V';
   std::complex<double>  worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;                         // for query of lwork
   double* rwork;
   int lrwork = 2 * Size;
   Fortran::integer info = 0;

   std::complex<double>* Data = new std::complex<double>[Size*Size];
   memcpy(Data, DataX, Size*Size*sizeof(std::complex<double>));

   rwork = new double[lrwork];

   // query for the optimial size of the workspace
   LAPACK::zgeev(jobl, jobr, Size, Data, LeadingDim, Eigen, LeftVectors, LeftLeadingDim,
                 RightVectors, RightLeadingDim, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zgeev(jobl, jobr, Size, Data, LeadingDim, Eigen, LeftVectors, LeftLeadingDim,
                 RightVectors, RightLeadingDim, work, lwork, rwork, info);

   CHECK(info == 0)("LAPACK::zgeev")(info);

   delete[] work;
   delete[] rwork;
   delete[] Data;
}

void GeneralizedEigenSymmetric(int Size, double* A, int ldA, double* B, int ldB,
                               int First, int Last, double* Eigenval, double* Z, int ldZ, double abstol)
{
   int itype = 1;  // mode 1: solve A*x = (lambda)*B*x
   char jobz = 'V';
   char range = 'I';
   char uplo = 'L';
   Fortran::integer M = Last-First+1;
   // Note: the LAPACK docs don't guarantee that only the first M elements of the eivenvalue array
   // are accessed, so to be safe we supply a new vector of length Size and do a copy
   double* W = new double[Size];
   double worksize;
   double* work = &worksize;
   int lwork = -1;                         // for query of lwork
   int liwork = 5 * Size;
   Fortran::integer* iwork;
   Fortran::integer* ifail;
   Fortran::integer info = 0;

   iwork = new Fortran::integer[liwork];
   ifail = new Fortran::integer[Size];

   // query for the optimial size of the workspace
   LAPACK::dsygvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
                  work, lwork, iwork, ifail, info);

   lwork = int(work[0]);
   work = new double[lwork];

   // do the actual call
   LAPACK::dsygvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
                  work, lwork, iwork, ifail, info);

   CHECK(info == 0)("LAPACK::dsygvx")(info);

   // copy the eigenvalues into the Eigenval array
   std::copy(W, W+M, Eigenval);

   delete[] work;
   delete[] ifail;
   delete[] iwork;
   delete[] W;
}

void GeneralizedEigenHermitian(int Size, std::complex<double>* A, int ldA, std::complex<double>* B, int ldB,
                               int First, int Last, double* Eigenval, std::complex<double>* Z, int ldZ, double abstol)
{
   int itype = 1;  // mode 1: solve A*x = (lambda)*B*x
   char jobz = 'V';
   char range = 'I';
   char uplo = 'L';
   Fortran::integer M = Last-First+1;
   // Note: the LAPACK docs don't guarantee that only the first M elements of the eivenvalue array
   // are accessed, so to be safe we supply a new vector of length Size and do a copy
   double* W = new double[Size];
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;                         // for query of lwork
   int lrwork = 7 * Size;
   double* rwork;
   int liwork = 5 * Size;
   Fortran::integer* iwork;
   Fortran::integer* ifail;
   Fortran::integer info = 0;

   rwork = new double[lrwork];
   iwork = new Fortran::integer[liwork];
   ifail = new Fortran::integer[Size];

   // query for the optimial size of the workspace
   LAPACK::zhegvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
                  work, lwork, rwork, iwork, ifail, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zhegvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
                  work, lwork, rwork, iwork, ifail, info);

   CHECK(info == 0)("LAPACK::zhegvx")(info);

   // copy the eigenvalues into the Eigenval array
   std::copy(W, W+M, Eigenval);

   delete[] work;
   delete[] ifail;
   delete[] iwork;
   delete[] rwork;
   delete[] W;
}

void SingularValueDecomposition(int Size1, int Size2, double* A, double* U,
                                double* D, double* VT)
{
   // For FORTRAN, we have to regard everything as the transpose.  The SDV of A^T is
   // VT^T * D * U^T, which means we have to swap all row/column numbers,
   // and also interchange VT and U.
   char jobu = 'S';
   char jobvt = 'S';
   int m = Size2;    // remembering that FORTRAN regards the matrices as the transpose, rows/cols are reversed.
   int n = Size1;
   int min_mn = std::min(m,n);
   double* a = A;
   int lda = Size2;
   double* s = D;
   double* u = VT;
   int ldu = Size2;
   double* vt = U;
   int ldvt = min_mn;
   double worksize;
   double* work = &worksize;
   int lwork = -1;
   Fortran::integer info = 0;

   // query for the optimial sizes of the workspace
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);

   lwork = int(work[0]);
   work = new double[lwork];

   // do the actual call
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
   CHECK(info == 0)("LAPACK::dgesvd")(info);

   delete[] work;
}

void SingularValueDecomposition(int Size1, int Size2,
                                std::complex<double>* A,
                                std::complex<double>* U,
                                double* D,
                                std::complex<double>* VH)
{
   // For FORTRAN, we have to regard everything as the transpose.  The SDV of A^T is
   // VH^T * D * U^T, which means we have to swap all row/column numbers,
   // and also interchange VH and U.
   int m = Size2;    // remembering that FORTRAN regards the matrices as the
                     // transpose, rows/cols are reversed.
   int n = Size1;
   int min_mn = std::min(m,n);
   std::complex<double>* a = A;
   int lda = Size2;
   double* s = D;
   std::complex<double>* u = VH;
   char jobu = u ? 'S' : 'N';
   int ldu = Size2;
   std::complex<double>* vh = U;
   char jobvh = vh ? 'S' : 'N';
   int ldvh = min_mn;
   DEBUG_CHECK(ldvh != 0); // This corner case is not allowed by LAPACK
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;
   int lrwork = 5 * min_mn;
   Fortran::integer info = 0;

   double* rwork = new double[lrwork];

   // query for the optimial sizes of the workspace
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);
   CHECK(info == 0)("LAPACK::zgesvd")(info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);
   CHECK(info == 0)("LAPACK::zgesvd")(info);

   delete[] work;
   delete[] rwork;
}

// let MAX = max(M,N)
// A is M x N
// U is M x MAX
// D is MAX x MAX
// VT is MAX x N
void SingularValueDecompositionFull(int Size1, int Size2, double* A, double* U,
                                    double* D, double* VT)
{
   // For FORTRAN, we have to regard everything as the transpose.  The SDV of A^T is
   // VT^T * D * U^T, which means we have to swap all row/column numbers,
   // and also interchange VT and U.
   int m = Size2;    // remembering that FORTRAN regards the matrices as the transpose, rows/cols are reversed.
   int n = Size1;
   int min_mn = std::min(m,n);
   int max_mn = std::max(m,n);
   double* a = A;
   int lda = m;
   double* s = D;
   double* u = VT;
   char jobu = u ? 'A' : 'N';
   int ldu = max_mn;
   double* vt = U;
   char jobvt = vt ? 'A' : 'N';
   int ldvt = n;
   double worksize;
   double* work = &worksize;
   int lwork = -1;
   Fortran::integer info = 0;

   // query for the optimial sizes of the workspace
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);

   lwork = int(work[0]);
   work = new double[lwork];

   // do the actual call
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
   CHECK(info == 0)("LAPACK::dgesvd")(info);

   // Only the first min_mn singular values are set; the remainder are zero
   std::memset(s + min_mn*sizeof(double), 0, (max_mn-min_mn)*sizeof(double));

   delete[] work;
}

// A is MxN
// U is MxM
// D is min(M,N)
// VH is NxN
//
// Recall Fortran uses column-major, so everything is transposed
// This means A^T = VH^T D U^T
// M -> n
// N -> m
// A'  = A^T is m x n
// U'  = VH^T is m x m
// D'  = D is max * max //min(m,n)
// VH' = U^T is n x n

void SingularValueDecompositionFull(int Size1, int Size2,
                                    std::complex<double>* A,
                                    std::complex<double>* U,
                                    double* D,
                                    std::complex<double>* VH)
{
   // For FORTRAN, we have to regard everything as the transpose.  The SDV of A^T is
   // VH^T * D * U^T, which means we have to swap all row/column numbers,
   // and also interchange VH and U.
   int m = Size2;    // remembering that FORTRAN regards the matrices as the transpose, rows/cols are reversed.
   int n = Size1;
   //int max_mn = std::max(m,n);
   int min_mn = std::min(m,n);
   int max_mn = std::max(m,n);
   std::complex<double>* a = A;
   int lda = m;
   double* s = D;
   std::complex<double>* u = VH;
   char jobu = u ? 'A' : 'N';
   int ldu = m;
   std::complex<double>* vh = U;
   char jobvh = vh ? 'A' : 'N';
   int ldvh = n;
   DEBUG_CHECK(ldvh != 0); // This corner case is not allowed by LAPACK
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;
   int lrwork = 5 * min_mn;
   Fortran::integer info = 0;

   double* rwork = new  double[lrwork];

   // query for the optimial sizes of the workspace
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);
   CHECK(info == 0)("LAPACK::zgesvd")(info);

   // Only the first min_mn singular values are set; the remainder are zero
   std::memset(s+min_mn, 0, (max_mn-min_mn)*sizeof(double));

   delete[] work;
   delete[] rwork;
}

void TridiagonalizeHermitian(int Size, std::complex<double>* A, int ldA,
                             double* Diag,
                             double* SubDiag)
{
   Fortran::integer info = 0;
   int lwork = -1;
   std::complex<double>* tau = new std::complex<double>[Size-1];
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;

   // query for the optimial sizes of the workspace
   LAPACK::zhetrd('U', Size, A, ldA, Diag, SubDiag, tau, work, lwork, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];

   // do the actual call
   LAPACK::zhetrd('U', Size, A, ldA, Diag, SubDiag, tau, work, lwork, info);
   CHECK(info == 0)("LAPACK::zhetrd")(info);

   delete[] work;
   delete[] tau;
}

void CholeskyUpper(int Size, std::complex<double>* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::zpotrf('U', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::zpotrf")(info);
}

void CholeskyUpper(int Size, double* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::dpotrf('U', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::zpotrf")(info);
}

void CholeskyLower(int Size, std::complex<double>* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::zpotrf('L', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::zpotrf")(info);
}

void CholeskyLower(int Size, double* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::dpotrf('L', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::dpotrf")(info);
}

void InvertHPD(int Size, std::complex<double>* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::zpotrf('U', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::zpotrf")(info);
   LAPACK::zpotri('U', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::zpotri")(info);

   // Copy the upper component into the lower
   for (int i = 0; i < Size; ++i)
   {
      for (int j = i+1; j < Size; ++j)
      {
         A[i*Size+j] = std::conj(A[j*Size+i]);
      }
   }
}

void InvertGeneral(int Size, std::complex<double>* A, int ldA)
{
   Fortran::integer info = 0;
   Fortran::integer* ipiv = new Fortran::integer[Size];
   LAPACK::zgetrf(Size , Size, A, ldA, ipiv, info);
   CHECK(info == 0)("LAPACK::zgetrf")(info);

   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;
   LAPACK::zgetri(Size, A, ldA, ipiv, work, lwork, info);

   lwork = int(work[0].real());
   work = new std::complex<double>[lwork];
   LAPACK::zgetri(Size, A, ldA, ipiv, work, lwork, info);
   CHECK(info == 0)("LAPACK::zgetri")(info);

   delete[] work;
   delete[] ipiv;
}

void InvertUpperTriangular(int Size, std::complex<double>* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::ztrtri('U', 'N', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::ztrtri")(info);
}

void InvertLowerTriangular(int Size, std::complex<double>* A, int ldA)
{
   Fortran::integer info = 0;
   LAPACK::ztrtri('L', 'N', Size, A, ldA, info);
   CHECK(info == 0)("LAPACK::ztrtri")(info);
}

void LQ_Factorize(int Size1, int Size2, double* A, int ldA, double* Tau)
{
   Fortran::integer info = 0;
   double worksize;
   double* Work = &worksize;
   int lWork = -1;
   LAPACK::dgelqf(Size1, Size2, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0]);
   Work = new double[lWork];
   LAPACK::dgelqf(Size1, Size2, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::dgelqf")(info);

   delete[] Work;
}

void LQ_Factorize(int Size1, int Size2, std::complex<double>* A, int ldA, std::complex<double>* Tau)
{
   Fortran::integer info = 0;
   std::complex<double> worksize;
   std::complex<double>* Work = &worksize;
   int lWork = -1;
   LAPACK::zgelqf(Size1, Size2, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0].real());
   Work = new std::complex<double>[lWork];
   LAPACK::zgelqf(Size1, Size2, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::zgelqf")(info);

   delete[] Work;
}

void LQ_Construct(int Size1, int Size2, int k, double* A, int ldA, double* Tau)
{
   DEBUG_CHECK(Size2 >= Size1 && Size1 >= 0)(Size1)(Size2);
   Fortran::integer info = 0;
   double worksize;
   double* Work = &worksize;
   int lWork = -1;
   LAPACK::dorglq(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0]);
   Work = new double[lWork];
   LAPACK::dorglq(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::dorglq")(info);

   delete[] Work;
}

void LQ_Construct(int Size1, int Size2, int k, std::complex<double>* A, int ldA, std::complex<double>* Tau)
{
   DEBUG_CHECK(Size2 >= Size1 && Size1 >= 0)(Size1)(Size2);
   Fortran::integer info = 0;
   std::complex<double> worksize;
   std::complex<double>* Work = &worksize;
   int lWork = -1;
   LAPACK::zunglq(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0].real());
   Work = new std::complex<double>[lWork];
   LAPACK::zunglq(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::zunglq")(info);

   delete[] Work;
}

void QR_Factorize(int Size1, int Size2, double* A, int ldA, double* Tau)
{
   Fortran::integer info = 0;
   double worksize;
   double* Work = &worksize;
   int lWork = -1;
   LAPACK::dgeqrf(Size1, Size2, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0]);
   Work = new double[lWork];
   LAPACK::dgeqrf(Size1, Size2, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::dgeqrf")(info);

   delete[] Work;
}

void QR_Factorize(int Size1, int Size2, std::complex<double>* A, int ldA, std::complex<double>* Tau)
{
   Fortran::integer info = 0;
   std::complex<double> worksize;
   std::complex<double>* Work = &worksize;
   int lWork = -1;
   LAPACK::zgeqrf(Size1, Size2, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0].real());
   Work = new std::complex<double>[lWork];
   LAPACK::zgeqrf(Size1, Size2, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::zgeqrf")(info);

   delete[] Work;
}

void QR_Construct(int Size1, int Size2, int k, double* A, int ldA, double* Tau)
{
   DEBUG_CHECK(Size1 >= Size2 && Size2 >= 0)(Size1)(Size2);
   Fortran::integer info = 0;
   double worksize;
   double* Work = &worksize;
   int lWork = -1;
   LAPACK::dorgqr(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0]);
   Work = new double[lWork];
   LAPACK::dorgqr(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::dorgqr")(info);

   delete[] Work;
}

void QR_Construct(int Size1, int Size2, int k, std::complex<double>* A, int ldA, std::complex<double>* Tau)
{
   DEBUG_CHECK(Size1 >= Size2 && Size2 >= 0)(Size1)(Size2);
   Fortran::integer info = 0;
   std::complex<double> worksize;
   std::complex<double>* Work = &worksize;
   int lWork = -1;
   LAPACK::zungqr(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0].real());
   Work = new std::complex<double>[lWork];
   LAPACK::zungqr(Size1, Size2, k, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::zungqr")(info);

   delete[] Work;
}


} // namespace Private

// QR

std::pair<Matrix<std::complex<double>>, Matrix<std::complex<double>>>
QR_Factorize(Matrix<std::complex<double>> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   CHECK(s1 >= s2)("For QR it is required that rows >= columns!")(s1)(s2);
   int sz = std::min(s1, s2);  // For the QR, we require that s1 >= s2, so we are guaranteed that sz==s2 here
   if (sz == 0)
   {
      return std::make_pair(Matrix<std::complex<double>>(s1,s2,0.0), Matrix<std::complex<double>>(s2,s2,0.0));
   }
   Vector<std::complex<double>> Tau(sz);
   Private::LQ_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*s2 matrix
   Matrix<std::complex<double>> Q(s1, s2, 0.0);
   Q(LinearAlgebra::all, LinearAlgebra::range(0,sz)) = M(LinearAlgebra::range(0,s1), LinearAlgebra::range(0,sz));

   Private::LQ_Construct(s2, s1, sz, data(Q), stride1(Q), data(Tau));

   // Zero the unused parts of m, which now becomes upper-triangular
   for (int i = 0; i < s2; ++i)
   {
      int msz = std::min(i,s2);
      for (int j = 0; j < msz; ++j)
      {
         M(i,j) = 0.0;
      }
   }
   return std::make_pair(std::move(Q), std::move(M));
}

std::pair<Matrix<double>, Matrix<double>>
QR_Factorize(Matrix<double> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   CHECK(s1 >= s2)("For QR it is required that rows >= columns!")(s1)(s2);
   int sz = std::min(s1, s2);
   if (sz == 0)
   {
      return std::make_pair(Matrix<double>(s1,s2,0.0), Matrix<double>(s2,s2,0.0));
   }
   Vector<double> Tau(sz);
   Private::LQ_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*s2 matrix
   Matrix<double> Q(s1, s2, 0.0);
   Q(LinearAlgebra::all, LinearAlgebra::range(0,sz)) = M(LinearAlgebra::range(0,s1), LinearAlgebra::range(0,sz));

   Private::LQ_Construct(s2, s1, sz, data(Q), stride1(Q), data(Tau));

   // Zero the unused parts of m, which now becomes upper-triangular
   for (int i = 0; i < s2; ++i)
   {
      int msz = std::min(i,s2);
      for (int j = 0; j < msz; ++j)
      {
         M(i,j) = 0.0;
      }
   }
   return std::make_pair(std::move(Q), std::move(M));
}
std::pair<Matrix<std::complex<double>>, Matrix<std::complex<double>>>
QR_FactorizeThin(Matrix<std::complex<double>> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   int sz = std::min(s1, s2);
   if (sz == 0)
   {
      return std::make_pair(Matrix<std::complex<double>>(s1,s2,0.0), Matrix<std::complex<double>>(s2,s2,0.0));
   }
   Vector<std::complex<double>> Tau(sz);
   Private::LQ_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*sz matrix
   Matrix<std::complex<double>> Q = M(LinearAlgebra::range(0,s1), LinearAlgebra::range(0,sz));

   Private::LQ_Construct(sz, s1, sz, data(Q), stride1(Q), data(Tau));

   // Copy the upper-triangular parts of M into R (new matrix, since it is a different size)
   Matrix<std::complex<double>> R(sz, s2, 0.0);
   for (int i = 0; i < sz; ++i)
   {
      for (int j = i; j < s2; ++j)
      {
         R(i,j) = M(i,j);
      }
   }
   return std::make_pair(std::move(Q), std::move(R));
}

std::pair<Matrix<double>, Matrix<double>>
QR_FactorizeThin(Matrix<double> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   int sz = std::min(s1, s2);
   if (sz == 0)
   {
      return std::make_pair(Matrix<double>(s1,sz,0.0), Matrix<double>(sz,s2,0.0));
   }
   Vector<double> Tau(sz);
   Private::LQ_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*sz matrix
   Matrix<double> Q = M(LinearAlgebra::range(0,s1), LinearAlgebra::range(0,sz));
   Private::LQ_Construct(sz, s1, sz, data(Q), stride1(Q), data(Tau));

   // Copy the upper-triangular parts of M into R (new matrix, since it is a different size)
   Matrix<double> R(sz, s2, 0.0);
   for (int i = 0; i < sz; ++i)
   {
      for (int j = i; j < s2; ++j)
      {
         R(i,j) = M(i,j);
      }
   }
   return std::make_pair(std::move(Q), std::move(R));
}

// LQ

std::pair<Matrix<std::complex<double>>, Matrix<std::complex<double>>>
LQ_Factorize(Matrix<std::complex<double>> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   int sz = std::min(s1, s2);  // For the LQ, we require that s2 >= s1, so we are guaranteed that sz==s1 here
   if (sz == 0)
   {
      return std::make_pair(Matrix<std::complex<double>>(s1,s2,0.0), Matrix<std::complex<double>>(s2,s2,0.0));
   }
   Vector<std::complex<double>> Tau(sz);
   Private::QR_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*s2 matrix
   Matrix<std::complex<double>> Q(s1, s2, 0.0);
   Q(LinearAlgebra::range(0,sz), LinearAlgebra::all) = M(LinearAlgebra::range(0,sz), LinearAlgebra::range(0,s2));
   Private::QR_Construct(s2, s1, sz, data(Q), stride1(Q), data(Tau));

   // Zero the unused parts of m, which now becomes lower-triangular
   for (int i = 0; i < sz; ++i)
   {
      for (int j = i+1; j < s2; ++j)
      {
         M(i,j) = 0.0;
      }
   }
   return std::make_pair(std::move(Q), std::move(M));
}

std::pair<Matrix<double>, Matrix<double>>
LQ_Factorize(Matrix<double> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   int sz = std::min(s1, s2);
   if (sz == 0)
   {
      return std::make_pair(Matrix<double>(s1,s2,0.0), Matrix<double>(s2,s2,0.0));
   }
   Vector<double> Tau(sz);
   Private::QR_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*s2 matrix
   Matrix<double> Q(s1, s2, 0.0);
   Q(LinearAlgebra::range(0,sz), LinearAlgebra::all) = M(LinearAlgebra::range(0,sz), LinearAlgebra::range(0,s2));
   Private::QR_Construct(s2, s1, sz, data(Q), stride1(Q), data(Tau));

   // Zero the unused parts of m, which now becomes lower-triangular
   for (int i = 0; i < sz; ++i)
   {
      for (int j = i+1; j < s2; ++j)
      {
         M(i,j) = 0.0;
      }
   }
   return std::make_pair(std::move(M), std::move(Q));
}
std::pair<Matrix<std::complex<double>>, Matrix<std::complex<double>>>
LQ_FactorizeThin(Matrix<std::complex<double>> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   int sz = std::min(s1, s2);
   if (sz == 0)
   {
      return std::make_pair(Matrix<std::complex<double>>(s1,s2,0.0), Matrix<std::complex<double>>(s2,s2,0.0));
   }
   Vector<std::complex<double>> Tau(sz);
   Private::QR_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*sz matrix
   Matrix<std::complex<double>> Q = M(LinearAlgebra::range(0,sz), LinearAlgebra::range(0,s2));
   Private::QR_Construct(s2, sz, sz, data(Q), stride1(Q), data(Tau));

   // Copy the lower-triangular parts of M into L (new matrix, since it is a different size)
   Matrix<std::complex<double>> L(s1, sz, 0.0);
   for (int i = 0; i < s1; ++i)
   {
      int msz = std::min(i+1,sz);
      for (int j = 0; j < msz; ++j)
      {
         L(i,j) = M(i,j);
      }
   }
   return std::make_pair(std::move(L), std::move(Q));
}

std::pair<Matrix<double>, Matrix<double>>
LQ_FactorizeThin(Matrix<double> M)
{
   int s1 = size1(M);
   int s2 = size2(M);
   int sz = std::min(s1, s2);
   if (sz == 0)
   {
      return std::make_pair(Matrix<double>(s1,sz,0.0), Matrix<double>(sz,s2,0.0));
   }
   Vector<double> Tau(sz);
   Private::QR_Factorize(size2(M), size1(M), data(M), stride1(M), data(Tau));

   // Convert the product of elementary reflectors into the Q matrix, as an s1*sz matrix
   Matrix<double> Q = M(LinearAlgebra::range(0,sz), LinearAlgebra::range(0,s2));
   Private::QR_Construct(s2, sz, sz, data(Q), stride1(Q), data(Tau));

   // Copy the lower-triangular parts of M into L (new matrix, since it is a different size)
   Matrix<double> L(s1, sz, 0.0);
   for (int i = 0; i < s1; ++i)
   {
      int msz = std::min(i+1,sz);
      for (int j = 0; j < msz; ++j)
      {
         L(i,j) = M(i,j);
      }
   }
   return std::make_pair(std::move(L), std::move(Q));
}

double amax(Matrix<double> const& X)
{
   double r = 0;
   for (int i = 0; i < X.size1(); ++i)
   {
      for (int j = 0; j < X.size2(); ++j)
      {
         if (std::abs(X(i,j)) < r || r == 0)
            r = std::abs(X(i,j));
      }
   }
   return r;
}

double amax(Matrix<std::complex<double>> const& X)
{
   double r2 = 0;
   for (int i = 0; i < X.size1(); ++i)
   {
      for (int j = 0; j < X.size2(); ++j)
      {
         auto x = std::norm(X(i,j));
         if (x < r2 || r2 == 0)
            r2 = x;
      }
   }
   return std::sqrt(r2);
}

Vector<std::complex<double>>
operator*(Matrix<std::complex<double>> const& M, Vector<std::complex<double>> const& v)
{
   Vector<std::complex<double>> Result(size1(M));
   for (unsigned i = 0; i < size1(M); ++i)
   {
      std::complex<double> x = 0.0;
      for (unsigned j = 0; j < size2(M); ++j)
      {
         x += M(i,j) * v[j];
      }
      Result[i] = x;
   }
   return Result;
}

Vector<std::complex<double>>
operator*(LinearAlgebra::MatrixTransposeProxy<LinearAlgebra::MatrixTransformProxy<const LinearAlgebra::Matrix<std::complex<double>, LinearAlgebra::RowMajor>&, LinearAlgebra::Conj<std::complex<double>>>> const& MM,
Vector<std::complex<double>> const& v)
{
   LinearAlgebra::Matrix<std::complex<double>> const& M = MM.base().base();
   Vector<std::complex<double>> Result(size2(M));
   for (unsigned i = 0; i < size2(M); ++i)
   {
      std::complex<double> x = 0.0;
      for (unsigned j = 0; j < size1(M); ++j)
      {
         x += conj(M(j,i)) * v[j];
      }
      Result[i] = x;
   }
   return Result;
}



std::pair<double, Vector<std::complex<double>>>
LeastSquaresRegularized(Matrix<std::complex<double>> const& A, Vector<std::complex<double>> const& b,
                        double alpha)
{
   Matrix<std::complex<double>> U, Vh;
   DiagonalMatrix<double> D;
   SingularValueDecompositionFullLeft(A, U, D, Vh);
   DiagonalMatrix<double> X = transform(D, [alpha](double x){ return x / (x*x + alpha);} );
   Vector<std::complex<double>> Result = herm(Vh) * (X * (herm(U) * b));
   double Residual = norm_2(U * ((D*X) * (herm(U) * Result)) - b);
   return std::make_pair(Residual, Result);
}

} // namespace LinearAlgebra
