// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/eigen.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// implementation of eigen.h LAPACK wrappers
//
// Created 2004-06-07 Ian McCulloch
//

#include "eigen.h"
#include "common/lapackf.h"
#include "common/stackallocator.h"
#include <algorithm>

namespace LinearAlgebra
{

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
   LWork = WorkSize.real();
   std::complex<double>* Work = new std::complex<double>[LWork];
   // do the actual call
   LAPACK::zgehrd(Size, ilo, ihi, Data, LeadingDim, Tau, Work, LWork, info);
   CHECK(info == 0)("LAPACK::zgehrd")(info);
   delete[] Work;

   // Third step: compute the eigenvalues of the Hessenberg matrix
   // workspace query
   std::complex<double> z;
   LAPACK::zhseqr('E', 'N', Size, ilo, ihi, Data, LeadingDim, Eigen, &z, 1, &WorkSize, -1, info);
   LWork = WorkSize.real();
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
   work = static_cast<double*>(StackAlloc::allocate(lwork * sizeof(double)));

   // do the actual call
   LAPACK::dsyev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, info);

   CHECK(info == 0)("LAPACK::dsyev")(info)(makeMatrix(Data, Size, LeadingDim));

   StackAlloc::deallocate(work, lwork * sizeof(double));
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

   rwork = static_cast<double*>(StackAlloc::allocate(lrwork * sizeof(double)));

   // query for the optimial size of the workspace
   LAPACK::zheev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = static_cast<std::complex<double>*>(StackAlloc::allocate(lwork * sizeof(std::complex<double>)));

   // do the actual call
   LAPACK::zheev(jobz, uplo, Size, Data, LeadingDim, Eigen, work, lwork, rwork, info);

   CHECK(info == 0)("LAPACK::zheev")(info);

   StackAlloc::deallocate(work , lwork * sizeof(std::complex<double>));
   StackAlloc::deallocate(rwork, lrwork * sizeof(double));
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

   std::complex<double>* Data = StackAlloc::allocate_type<std::complex<double> >(Size*Size);
   memcpy(Data, DataX, Size*Size*sizeof(std::complex<double>));

   rwork = static_cast<double*>(StackAlloc::allocate(lrwork * sizeof(double)));

   // query for the optimial size of the workspace
   LAPACK::zgeev(jobl, jobr, Size, Data, LeadingDim, Eigen, LeftVectors, LeftLeadingDim, 
                 RightVectors, RightLeadingDim, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = static_cast<std::complex<double>*>(StackAlloc::allocate(lwork * sizeof(std::complex<double>)));

   // do the actual call
   LAPACK::zgeev(jobl, jobr, Size, Data, LeadingDim, Eigen, LeftVectors, LeftLeadingDim, 
                 RightVectors, RightLeadingDim, work, lwork, rwork, info);

   CHECK(info == 0)("LAPACK::zgeev")(info);

   StackAlloc::deallocate(work , lwork * sizeof(std::complex<double>));
   StackAlloc::deallocate(rwork, lrwork * sizeof(double));
   StackAlloc::deallocate_type<std::complex<double> >(Data, Size*Size);
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
   double* W = static_cast<double*>(StackAlloc::allocate(Size * sizeof(double)));
   double worksize;
   double* work = &worksize;
   int lwork = -1;                         // for query of lwork
   int liwork = 5 * Size;
   Fortran::integer* iwork;
   Fortran::integer* ifail;
   Fortran::integer info = 0;

   iwork = static_cast<Fortran::integer*>(StackAlloc::allocate(liwork * sizeof(Fortran::integer)));
   ifail = static_cast<Fortran::integer*>(StackAlloc::allocate(Size * sizeof(Fortran::integer)));

   // query for the optimial size of the workspace
   LAPACK::dsygvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
		  work, lwork, iwork, ifail, info);

   lwork = int(work[0]);
   work = static_cast<double*>(StackAlloc::allocate(lwork * sizeof(double)));

   // do the actual call
   LAPACK::dsygvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
		  work, lwork, iwork, ifail, info);

   CHECK(info == 0)("LAPACK::dsygvx")(info);

   // copy the eigenvalues into the Eigenval array
   std::copy(W, W+M, Eigenval);

   StackAlloc::deallocate(work, lwork * sizeof(double));
   StackAlloc::deallocate(ifail, Size * sizeof(Fortran::integer));
   StackAlloc::deallocate(iwork, liwork * sizeof(Fortran::integer));
   StackAlloc::deallocate(W, Size * sizeof(double));
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
   double* W = static_cast<double*>(StackAlloc::allocate(Size * sizeof(double)));
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;                         // for query of lwork
   int lrwork = 7 * Size;
   double* rwork;
   int liwork = 5 * Size;
   Fortran::integer* iwork;
   Fortran::integer* ifail;
   Fortran::integer info = 0;

   rwork = static_cast<double*>(StackAlloc::allocate(lrwork * sizeof(double)));
   iwork = static_cast<Fortran::integer*>(StackAlloc::allocate(liwork * sizeof(Fortran::integer)));
   ifail = static_cast<Fortran::integer*>(StackAlloc::allocate(Size * sizeof(Fortran::integer)));

   // query for the optimial size of the workspace
   LAPACK::zhegvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
		  work, lwork, rwork, iwork, ifail, info);

   lwork = int(work[0].real());
   work = static_cast<std::complex<double>*>(StackAlloc::allocate(lwork * sizeof(std::complex<double>)));

   // do the actual call
   LAPACK::zhegvx(itype, jobz, range, uplo, Size, A, ldA, B, ldB, 0, 0, First, Last, abstol, M, W, Z, ldZ,
		  work, lwork, rwork, iwork, ifail, info);

   CHECK(info == 0)("LAPACK::zhegvx")(info);

   // copy the eigenvalues into the Eigenval array
   std::copy(W, W+M, Eigenval);

   StackAlloc::deallocate(work, lwork * sizeof(std::complex<double>));
   StackAlloc::deallocate(ifail, Size * sizeof(Fortran::integer));
   StackAlloc::deallocate(iwork, liwork * sizeof(Fortran::integer));
   StackAlloc::deallocate(rwork, lrwork * sizeof(double));
   StackAlloc::deallocate(W, Size * sizeof(double));
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
   work = static_cast<double*>(StackAlloc::allocate(lwork * sizeof(double)));

   // do the actual call
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
   CHECK(info == 0)("LAPACK::dgesvd")(info);

   StackAlloc::deallocate(work, lwork * sizeof(double));
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
   char jobu = 'S';
   char jobvh = 'S';
   int m = Size2;    // remembering that FORTRAN regards the matrices as the 
                     // transpose, rows/cols are reversed.
   int n = Size1;
   int min_mn = std::min(m,n);
   std::complex<double>* a = A;
   int lda = Size2;
   double* s = D;
   std::complex<double>* u = VH;
   int ldu = Size2;
   std::complex<double>* vh = U;
   int ldvh = min_mn;
   DEBUG_CHECK(ldvh != 0); // This corner case is not allowed by LAPACK
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;
   int lrwork = 5 * min_mn;
   Fortran::integer info = 0;

   double* rwork = StackAlloc::allocate_type<double>(lrwork);

   // query for the optimial sizes of the workspace
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = StackAlloc::allocate_type<std::complex<double> >(lwork);

   // do the actual call
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);
   CHECK(info == 0)("LAPACK::zgesvd")(info);

   StackAlloc::deallocate_type<std::complex<double> >(work, lwork);
   StackAlloc::deallocate_type<double>(rwork, lrwork);
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
   char jobu = 'A';
   char jobvt = 'A';
   int m = Size2;    // remembering that FORTRAN regards the matrices as the transpose, rows/cols are reversed.
   int n = Size1;
   int max_mn = std::max(m,n);
   double* a = A;
   int lda = m;
   double* s = D;
   double* u = VT;
   int ldu = max_mn;
   double* vt = U;
   int ldvt = n;
   double worksize;
   double* work = &worksize;
   int lwork = -1;
   Fortran::integer info = 0;

   // query for the optimial sizes of the workspace
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);

   lwork = int(work[0]);
   work = static_cast<double*>(StackAlloc::allocate(lwork * sizeof(double)));

   // do the actual call
   LAPACK::dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
   CHECK(info == 0)("LAPACK::dgesvd")(info);

   StackAlloc::deallocate(work, lwork * sizeof(double));
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
// D'  = D is min(m,n)
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
   char jobu = 'A';
   char jobvh = 'A';
   int m = Size2;    // remembering that FORTRAN regards the matrices as the transpose, rows/cols are reversed.
   int n = Size1;
   int max_mn = std::max(m,n);
   int min_mn = std::min(m,n);
   std::complex<double>* a = A;
   int lda = m;
   double* s = D;
   std::complex<double>* u = VH;
   int ldu = m;
   std::complex<double>* vh = U;
   int ldvh = n;
   DEBUG_CHECK(ldvh != 0); // This corner case is not allowed by LAPACK
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;
   int lrwork = 5 * min_mn;
   Fortran::integer info = 0;

   double* rwork = StackAlloc::allocate_type<double>(lrwork);

   // query for the optimial sizes of the workspace
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);

   lwork = int(work[0].real());
   work = StackAlloc::allocate_type<std::complex<double> >(lwork);

   // do the actual call
   LAPACK::zgesvd(jobu, jobvh, m, n, a, lda, s, u, ldu, vh, ldvh, work, lwork, rwork, info);
   CHECK(info == 0)("LAPACK::zgesvd")(info);

   StackAlloc::deallocate_type<std::complex<double> >(work, lwork);
   StackAlloc::deallocate_type<double>(rwork, lrwork);
}

void TridiagonalizeHermitian(int Size, std::complex<double>* A, int ldA, 
                             double* Diag,
                             double* SubDiag)
{
   Fortran::integer info = 0;
   int lwork = -1;
   std::complex<double>* tau = StackAlloc::allocate_type<std::complex<double> >(Size-1);
   std::complex<double> worksize;
   std::complex<double>* work = &worksize;

   // query for the optimial sizes of the workspace
   LAPACK::zhetrd('U', Size, A, ldA, Diag, SubDiag, tau, work, lwork, info);

   lwork = int(work[0].real());
   work = StackAlloc::allocate_type<std::complex<double> >(lwork);

   // do the actual call
   LAPACK::zhetrd('U', Size, A, ldA, Diag, SubDiag, tau, work, lwork, info);
   CHECK(info == 0)("LAPACK::zhetrd")(info);

   StackAlloc::deallocate_type<std::complex<double> >(work, lwork);
   StackAlloc::deallocate_type<std::complex<double> >(tau, Size-1);
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
   Fortran::integer* ipiv = StackAlloc::allocate_type<Fortran::integer>(Size);
   LAPACK::zgetrf(Size , Size, A, ldA, ipiv, info);
   CHECK(info == 0)("LAPACK::zgetrf")(info);

   std::complex<double> worksize;
   std::complex<double>* work = &worksize;
   int lwork = -1;
   LAPACK::zgetri(Size, A, ldA, ipiv, work, lwork, info);

   lwork = int(work[0].real());
   work = StackAlloc::allocate_type<std::complex<double> >(lwork);
   LAPACK::zgetri(Size, A, ldA, ipiv, work, lwork, info);
   CHECK(info == 0)("LAPACK::zgetri")(info);

   StackAlloc::deallocate_type<std::complex<double> >(work, lwork);
   StackAlloc::deallocate_type<Fortran::integer>(ipiv, Size);
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

void LQ_Factorize(int Size1, int Size2, std::complex<double>* A, int ldA, std::complex<double>* Tau)
{
   Fortran::integer info = 0;
   std::complex<double> worksize;
   std::complex<double>* Work = &worksize;
   int lWork = -1;
   LAPACK::zgelqf(Size1, Size2, A, ldA, Tau, Work, lWork, info);

   lWork = int(Work[0].real());
   Work = StackAlloc::allocate_type<std::complex<double> >(lWork);
   LAPACK::zgelqf(Size1, Size2, A, ldA, Tau, Work, lWork, info);
   CHECK(info == 0)("LAPACK::zgelqf")(info);

   StackAlloc::deallocate_type<std::complex<double> >(Work, lWork);
}

} // namespace Private

} // namespace LinearAlgebra
