// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/dataops_product.h
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
//
// dataops_product.h
//
// primitive operations for calculating products
//

#if !defined(DATAOPS_PRODUCT_H_HCV3478YT78YGHLUIH3Q789YJFLORJET9O)
#define DATAOPS_PRODUCT_H_HCV3478YT78YGHLUIH3Q789YJFLORJET9O

namespace ops
{

// in the absence of a decent multiply strategy, we have the following bad, tempoary generating functions.
// multiply seems to be hard to make fast in the general case - maybe we should just make a heap of these
// sorts of functions and patch them all together into an ET heirachy?
template <class Scalar, class D1, class D2>
DenseMatrix<Scalar> product(RandomAccessMatrix<Scalar const, D1> const& M1, RandomAccessMatrix<Scalar const, D2> const& M2)
{
   PRECONDITION(M1.cols() == M2.rows());
   int Rows = M1.rows();
   int Cols = M2.cols();
   int Inner = M1.cols();
   DenseMatrix<Scalar> Result(Rows, Cols);
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(i,k) * M2(k,j);
         Result(i,j) = Acc;
      }
   }
   return Result;
}

#if 0

// product A * transpose(B)
template <class Scalar>
DenseMatrix<Scalar> product_ABt(DenseMatrix<Scalar> const& M1, DenseMatrix<Scalar> const& M2)
{
   PRECONDITION(M1.cols() == M2.cols());
   int Rows = M1.rows();
   int Cols = M2.rows();
   int Inner = M1.cols();
   DenseMatrix<Scalar> Result(Rows, Cols);
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(i,k) * M2(j,k);
         Result(i,j) = Acc;
      }
   }
   return Result;
}

template <class Scalar>
DenseMatrix<Scalar> product_AtB(DenseMatrix<Scalar> const& M1, DenseMatrix<Scalar> const& M2)
{
   PRECONDITION(M1.rows() == M2.rows());
   int Rows = M1.cols();
   int Cols = M2.cols();
   int Inner = M1.rows();
   DenseMatrix<Scalar> Result(Rows, Cols);
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(k,i) * M2(k,j);
         Result(i,j) = Acc;
      }
   }
   return Result;
}

#endif

// add product A * B
template <class Scalar, class D2, class D3>
void add_product_AB(DenseMatrix<Scalar>& Result,
                    RandomAccessMatrix<Scalar const, D2> const& M1,
                    RandomAccessMatrix<Scalar const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.rows());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.cols());
   int Rows = M1.rows();
   int Cols = M2.cols();
   int Inner = M1.cols();
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(i,k) * M2(k,j);
         Result(i,j) += Acc;
      }
   }
}

template <class Scalar, class D2, class D3>
void add_product_sAB(DenseMatrix<Scalar>& Result,
                     Scalar s,
                     RandomAccessMatrix<Scalar const, D2> const& M1,
                     RandomAccessMatrix<Scalar const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.rows());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.cols());
   int Rows = M1.rows();
   int Cols = M2.cols();
   int Inner = M1.cols();
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(i,k) * M2(k,j);
         Result(i,j) += s * Acc;
      }
   }
}

template <class Scalar, class D2, class D3>
void add_product_ABt(DenseMatrix<Scalar>& Result,
                     RandomAccessMatrix<Scalar const, D2> const& M1,
                     RandomAccessMatrix<Scalar const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.cols());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.rows());
   int Rows = M1.rows();
   int Cols = M2.rows();
   int Inner = M1.cols();
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(i,k) * M2(j,k);
         Result(i,j) += Acc;
      }
   }
}

// add product s * A * transpose(B)
template <class Scalar, class D2, class D3>
void add_product_sABt(DenseMatrix<Scalar>& Result,
                      Scalar s,
                      RandomAccessMatrix<Scalar const, D2> const& M1,
                      RandomAccessMatrix<Scalar const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.cols());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.rows());
   int Rows = M1.rows();
   int Cols = M2.rows();
   int Inner = M1.cols();
   for (int i = 0; i < Rows; ++i)
   {
      for (int j = 0; j < Cols; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < Inner; ++k)
            Acc += M1(i,k) * M2(j,k);
         Result(i,j) += s * Acc;
      }
   }
}

#if 0

// add product of s * A * B * transpose(C) to Result
template <class Scalar, class D2, class D3, class D4>
void add_product_sABCt(DenseMatrix<Scalar>& Result,
                       Scalar s,
                       RandomAccessMatrix<Scalar const, D2> const& A,
                       RandomAccessMatrix<Scalar const, D3> const& B,
                       RandomAccessMatrix<Scalar const, D4> const& C)
{
   // implements R(i,l) = A(i,j) * B(j,k) * C(l,k)
   PRECONDITION(Result.rows() == A.rows());
   PRECONDITION(A.cols() == B.rows());
   PRECONDITION(B.cols() == C.cols());
   PRECONDITION(Result.cols() == C.rows());
   int const isz = Result.rows();
   int const lsz = Result.cols();
   int const jsz = A.cols();
   int const ksz = B.cols();
   Vector<Scalar> TempBC(jsz);

   for (int l = 0; l < lsz; ++l)
   {
      for (int j = 0; j < jsz; ++j)
      {
         Scalar Acc = 0;
         for (int k = 0; k < ksz; ++k)
         {
            Acc += B(j,k) * C(l,k);
         }
         TempBC[j] = Acc;
      }
      for (int i = 0; i < isz; ++i)
      {
         Scalar Acc = 0;
         for (int j = 0; j < jsz; ++j)
         {
            Acc += A(i,j) * TempBC[j];
         }
         Result(i,l) += s * Acc;
      }
   }
}

#endif

#include <iostream>

// add product of s * A * B * transpose(C) to Result
template <class Scalar, class D2, class D3, class D4>
void add_product_sABCt(DenseMatrix<Scalar>& Result,
                       Scalar s,
                       Stride1Matrix<Scalar const, D2> const& A,
                       Stride1Matrix<Scalar const, D3> const& B,
                       Stride1Matrix<Scalar const, D4> const& C)
{
   // implements R(i,l) = A(i,j) * B(j,k) * C(l,k)
   PRECONDITION(Result.rows() == A.rows());
   PRECONDITION(A.cols() == B.rows());
   PRECONDITION(B.cols() == C.cols());
   PRECONDITION(Result.cols() == C.rows());
   int const isz = Result.rows();
   int const lsz = Result.cols();
   int const jsz = A.cols();
   int const ksz = B.cols();
   int const ldR = Result.leading_dimension();
   int const ldA = A.leading_dimension();
   int const ldB = B.leading_dimension();
   int const ldC = C.leading_dimension();
   Vector<Scalar> TempBC(jsz);

   Scalar* ResultPtr = Result.data();
   Scalar const* CPtr = C.data();
   for (int l = 0; l < lsz; ++l)
   {
      Scalar* TP = &TempBC[0];
      Scalar const* BPtr = B.data();
      for (int j = 0; j < jsz; ++j)
      {
         Scalar const* BP = BPtr;
         Scalar Acc = 0;
         Scalar const* CP = CPtr;
         for (int k = 0; k < ksz; ++k)
         {
            Acc += (*BP++) * (*CP++);  // Acc += B(j,k) * C(l,k)
         }
         (*TP++) = Acc;      // TempBC[j] = Acc
         BPtr += ldB;        // next column of B
      }

      Scalar* RP = ResultPtr;
      Scalar const* APtr = A.data();
      for (int i = 0; i < isz; ++i)
      {
         Scalar const* AP = APtr;
         Scalar const* TP = &TempBC[0];
         Scalar Acc = 0;
         for (int j = 0; j < jsz; ++j)
         {
            Acc += (*AP++) * (*TP++);   // Acc += A(i,j) * TempBC[j]
         }
         *RP += s * Acc;   // Result(i,l) += s * Acc
         RP += ldR;        // next row of R
         APtr += ldA;
      }
      CPtr += ldC;
      ++ResultPtr; // next column of R
   }
}

#include "common/blas3f.h"

inline
void product_sABt_blas(int isz, int jsz, int ksz, double* restrict R, double s, double const* A, double const* B)
{
   BLAS::dgemm('T', 'N', jsz, isz, ksz, s, B, ksz, A, ksz, 0, R, jsz);
}

inline
void product_sAB_blas(int isz, int jsz, int ksz, double* restrict R, double s, double const* A, double const* B)
{
   BLAS::dgemm('N', 'N', jsz, isz, ksz, s, B, jsz, A, ksz, 0, R, jsz);
}

inline
void product_sAtB_blas(int A_cols, int B_rows, int B_cols, double* restrict R, double s, double const* A, double const* B)
{
   BLAS::dgemm('N', 'T', B_cols, A_cols, B_rows, s, B, B_cols, A, A_cols, 0, R, B_cols);
}

inline
void product_AB_blas(int isz, int jsz, int ksz, double* restrict R, double const* A, double const* B)
{
   BLAS::dgemm('N', 'N', jsz, isz, ksz, 1.0, B, jsz, A, ksz, 0, R, jsz);
}

inline
void product_ABt_blas(int isz, int jsz, int ksz, double* restrict R, double const* A, double const* B)
{
   BLAS::dgemm('T', 'N', jsz, isz, ksz, 1.0, B, ksz, A, ksz, 0, R, jsz);
}

inline
void product_AtB_blas(int A_cols, int B_rows, int B_cols, double* restrict R, double const* A, double const* B)
{
   BLAS::dgemm('N', 'T', B_cols, A_cols, B_rows, 1.0, B, B_cols, A, A_cols, 0, R, B_cols);
}

inline
void add_product_AB_blas(int isz, int jsz, int ksz, double* restrict R, double const* A, double const* B)
{
   BLAS::dgemm('N', 'N', jsz, isz, ksz, 1.0, B, jsz, A, ksz, 1.0, R, jsz);
}

inline
void add_product_ABt_blas(int isz, int jsz, int ksz, double* restrict R, double const* A, double const* B)
{
   BLAS::dgemm('T', 'N', jsz, isz, ksz, 1.0, B, ksz, A, ksz, 1.0, R, jsz);
}

inline
void add_product_sAB_blas(int isz, int jsz, int ksz, double* restrict R, double s, double const* A, double const* B)
{
   BLAS::dgemm('N', 'N', jsz, isz, ksz, s, B, jsz, A, ksz, 1.0, R, jsz);
}

inline
void add_product_sABt_blas(int isz, int jsz, int ksz, double* restrict R, double s, double const* A, double const* B)
{
   BLAS::dgemm('T', 'N', jsz, isz, ksz, s, B, ksz, A, ksz, 1.0, R, jsz);
}

inline
void add_product_sAtB_blas(int A_cols, int B_rows, int B_cols, double* restrict R, double s, double const* A, double const* B)
{
   BLAS::dgemm('N', 'T', B_cols, A_cols, B_rows, s, B, B_cols, A, A_cols, 1.0, R, B_cols);
}

template <class D2, class D3, class D4>
void product_sABCt(DenseMatrix<double>& Result,
                       double s,
                       Stride1Matrix<double const, D2> const& A,
                       Stride1Matrix<double const, D3> const& B,
                       Stride1Matrix<double const, D4> const& C)
{
  // does R(i,l) = s * A(i,j) * B(j,k) * C(l,k)
   PRECONDITION(Result.rows() == A.rows());
   PRECONDITION(A.cols() == B.rows());
   PRECONDITION(B.cols() == C.cols());
   PRECONDITION(Result.cols() == C.rows());
   int const isz = Result.rows();
   int const lsz = Result.cols();
   int const jsz = A.cols();
   int const ksz = B.cols();
   int const ldR = Result.leading_dimension();
   int const ldA = A.leading_dimension();
   int const ldB = B.leading_dimension();
   int const ldC = C.leading_dimension();

   // order the multiplies so we do the minimum operations
   if ((ksz+isz)*(jsz*ksz) <= (ksz*isz)*(jsz+ksz))
   {
      Vector<double> TempBC(jsz * lsz);
      product_sABt_blas(jsz, lsz, ksz, &TempBC[0], s, B.data(), C.data());
      product_AB_blas(isz, lsz, jsz, Result.data(), A.data(), &TempBC[0]);
   }
   else
   {
      Vector<double> TempBC(isz * ksz);
      product_sAB_blas(isz, ksz, jsz, &TempBC[0], s, A.data(), B.data());
      product_ABt_blas(isz, lsz, ksz, Result.data(), &TempBC[0], C.data());
   }
}

template <class D2, class D3, class D4>
void add_product_sABCt(DenseMatrix<double>& Result,
                       double s,
                       Stride1Matrix<double const, D2> const& A,
                       Stride1Matrix<double const, D3> const& B,
                       Stride1Matrix<double const, D4> const& C)
{
  // does R(i,l) = s * A(i,j) * B(j,k) * C(l,k)
   PRECONDITION(Result.rows() == A.rows());
   PRECONDITION(A.cols() == B.rows());
   PRECONDITION(B.cols() == C.cols());
   PRECONDITION(Result.cols() == C.rows());
   int const isz = Result.rows();
   int const lsz = Result.cols();
   int const jsz = A.cols();
   int const ksz = B.cols();
   int const ldR = Result.leading_dimension();
   int const ldA = A.leading_dimension();
   int const ldB = B.leading_dimension();
   int const ldC = C.leading_dimension();

   // order the multiplies so we do the minimum operations
   if ((ksz+isz)*(jsz*ksz) <= (ksz*isz)*(jsz+ksz))
   {
      Vector<double> TempBC(jsz * lsz);
      product_sABt_blas(jsz, lsz, ksz, &TempBC[0], s, B.data(), C.data());
      add_product_AB_blas(isz, lsz, jsz, Result.data(), A.data(), &TempBC[0]);
   }
   else
   {
      Vector<double> TempBC(isz * ksz);
      product_sAB_blas(isz, ksz, jsz, &TempBC[0], s, A.data(), B.data());
      add_product_ABt_blas(isz, lsz, ksz, Result.data(), &TempBC[0], C.data());
   }
}

// This version uses a caller supplied buffer to store the tempoary results.
// The buffer size must be at least min(B.rows() * C.rows(), A.rows() * B.cols())
template <class D2, class D3, class D4>
void add_product_sABCt(DenseMatrix<double>& Result,
                       double s,
                       Stride1Matrix<double const, D2> const& A,
                       Stride1Matrix<double const, D3> const& B,
                       Stride1Matrix<double const, D4> const& C,
                       double* TempBC)
{
  // does R(i,l) = s * A(i,j) * B(j,k) * C(l,k)
   PRECONDITION(Result.rows() == A.rows());
   PRECONDITION(A.cols() == B.rows());
   PRECONDITION(B.cols() == C.cols());
   PRECONDITION(Result.cols() == C.rows());
   int const isz = Result.rows();
   int const lsz = Result.cols();
   int const jsz = A.cols();
   int const ksz = B.cols();
   int const ldR = Result.leading_dimension();
   int const ldA = A.leading_dimension();
   int const ldB = B.leading_dimension();
   int const ldC = C.leading_dimension();

#if 0

   typedef double Scalar;

   Scalar* ResultPtr = Result.data();
   Scalar const* CPtr = C.data();
   for (int l = 0; l < lsz; ++l)
   {
      Scalar* TP = &TempBC[0];
      Scalar const* BPtr = B.data();
      for (int j = 0; j < jsz; ++j)
      {
         Scalar const* BP = BPtr;
         Scalar Acc = 0;
         Scalar const* CP = CPtr;
         for (int k = 0; k < ksz; ++k)
         {
            Acc += (*BP++) * (*CP++);  // Acc += B(j,k) * C(l,k)
         }
         (*TP++) = Acc;      // TempBC[j] = Acc
         BPtr += ldB;        // next column of B
      }

      Scalar* RP = ResultPtr;
      Scalar const* APtr = A.data();
      for (int i = 0; i < isz; ++i)
      {
         Scalar const* AP = APtr;
         Scalar const* TP = &TempBC[0];
         Scalar Acc = 0;
         for (int j = 0; j < jsz; ++j)
         {
            Acc += (*AP++) * (*TP++);   // Acc += A(i,j) * TempBC[j]
         }
         *RP += s * Acc;   // Result(i,l) += s * Acc
         RP += ldR;        // next row of R
         APtr += ldA;
      }
      CPtr += ldC;
      ++ResultPtr; // next column of R
   }

#else

   // order the multiplies so we do the minimum operations
   if ((ksz+isz)*(jsz*ksz) <= (ksz*isz)*(jsz+ksz))
   {
      product_sABt_blas(jsz, lsz, ksz, &TempBC[0], s, B.data(), C.data());
      add_product_AB_blas(isz, lsz, jsz, Result.data(), A.data(), &TempBC[0]);
   }
   else
   {
      product_sAB_blas(isz, ksz, jsz, &TempBC[0], s, A.data(), B.data());
      add_product_ABt_blas(isz, lsz, ksz, Result.data(), &TempBC[0], C.data());
   }

#endif
}

template <class D2, class D3>
void add_product_AB(DenseMatrix<double>& Result,
                    Stride1Matrix<double const, D2> const& M1,
                    Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.rows());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.cols());
   add_product_AB_blas(M1.rows(), M2.cols(), M1.cols(), Result.data(), M1.data(), M2.data());
}

template <class D2, class D3>
void add_product_sAB(DenseMatrix<double>& Result,
                     double s,
                     Stride1Matrix<double const, D2> const& M1,
                     Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.rows());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.cols());
   add_product_sAB_blas(M1.rows(), M2.cols(), M1.cols(), Result.data(), s, M1.data(), M2.data());
}

template <class D2, class D3>
void add_product_ABt(DenseMatrix<double>& Result,
                     Stride1Matrix<double const, D2> const& M1,
                     Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.cols());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.rows());
   add_product_ABt_blas(M1.rows(), M2.rows(), M1.cols(), Result.data(), M1.data(), M2.data());
}

template <class D2, class D3>
void add_product_sABt(DenseMatrix<double>& Result,
                      double s,
                      Stride1Matrix<double const, D2> const& M1,
                      Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.cols());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.rows());
   add_product_sABt_blas(M1.rows(), M2.rows(), M1.cols(), Result.data(), s, M1.data(), M2.data());
}

template <class D2, class D3>
void add_product_sAtB(DenseMatrix<double>& Result,
                      double s,
                      Stride1Matrix<double const, D2> const& M1,
                      Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.rows() == M2.rows());
   PRECONDITION(Result.rows() == M1.cols());
   PRECONDITION(Result.cols() == M2.cols());
   add_product_sAtB_blas(M1.cols(), M2.rows(), M2.cols(), Result.data(), s, M1.data(), M2.data());
}

template <class D2, class D3>
void product_sABt(DenseMatrix<double>& Result,
                  double s,
                  Stride1Matrix<double const, D2> const& M1,
                  Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.cols() == M2.cols());
   PRECONDITION(Result.rows() == M1.rows());
   PRECONDITION(Result.cols() == M2.rows());
   product_sABt_blas(M1.rows(), M2.rows(), M1.cols(), Result.data(), s, M1.data(), M2.data());
}

template <class D2, class D3>
void product_sAtB(DenseMatrix<double>& Result,
                      double s,
                      Stride1Matrix<double const, D2> const& M1,
                      Stride1Matrix<double const, D3> const& M2)
{
   PRECONDITION(M1.rows() == M2.rows());
   PRECONDITION(Result.rows() == M1.cols());
   PRECONDITION(Result.cols() == M2.cols());
   product_sAtB_blas(M1.cols(), M2.rows(), M2.cols(), Result.data(), s, M1.data(), M2.data());
}

} // namespace ops

#endif
