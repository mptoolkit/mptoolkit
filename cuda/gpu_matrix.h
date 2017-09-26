// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cublas/gpu_matrix.h
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

#if !defined(MPTOOLKIT_CUBLAS_GPU_MATRIX_H)
#define MPTOOLKIT_CUBLAS_GPU_MATRIX_H

#include "cublas.h"
#include "gpu_buffer.h"
#include "blas/matrixref.h"
#include "blas/matrix.h"

namespace cublas
{

//
// GPU-storage matrix type.  Column-major format for compatability with BLAS
//

template <typename T>
class gpu_matrix : public blas::BlasMatrix<T, gpu_matrix<T>>
{
   public:
      typedef T value_type;

      gpu_matrix(int Rows_, int Cols_, arena const& A);

      gpu_matrix(gpu_matrix&& Other) = default;

      gpu_matrix(gpu_matrix const&) = delete;

      gpu_matrix& operator=(gpu_matrix&&) = delete;

      ~gpu_matrix()
      {
      }

      gpu_matrix& operator=(blas::Matrix<T> const& Other)
      {
         assign(*this, Other);
      }

      gpu_matrix& operator=(gpu_matrix const& Other)
      {
         assign(*this, Other);
      }

      // assignment of expressions based on the same matrix type
      template <typename U>
      gpu_matrix& operator=(blas::MatrixRef<T, gpu_matrix<T>, U> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      constexpr char trans() const { return 'N'; }

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dimension() const { return LeadingDimension; }

      cuda::gpu_buffer<T>& buffer() { return Buf; }
      cuda::gpu_buffer<T> const& buffer() const { return Buf; }

   private:
      int Rows;
      int Cols;
      int LeadingDimension;
      cuda::gpu_buffer<T> Buf;
};

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, arena const& Arena)
   : Rows(Rows_), Cols(Cols_), LeadingDimension(Rows_),
     Buf(cuda::gpu_buffer<T>::allocate(Rows_*Cols_, Arena))
{
}

// blocking matrix get
template <typename T>
blas::Matrix<T>
get_wait(gpu_matrix<T> const& M)
{
   blas::Matrix<T> Result(M.rows(), M.cols());
   cublas::check_error(cublasGetMatrix(M.rows(), M.cols(), sizeof(T),
				       M.buffer().device_ptr(), M.leading_dimension(),
				       Result.data(), Result.leading_dimension()));
   return Result;
}

template <typename T>
void
set_wait(gpu_matrix<T>& A, blas::Matrix<T> const& B)
{
   //TRACE(A.cols())(A.rows())(A.device_ptr())(A.leading_dim())(B.data())(leading_dimension(B));
   cublas::check_error(cublasSetMatrix(A.rows(), A.cols(), sizeof(T),
				       B.data(), leading_dimension(B),
				       A.device_ptr(), A.leading_dim()));
}

template <typename T>
void
assign(gpu_matrix<T>& A, blas::Matrix<T> const& B)
{
   cublas::check_error(cublasSetMatrixAsync(A.rows(), A.cols(), sizeof(T),
                                            B.data(), B.leading_dimension(),
                                            A.buffer().device_ptr(), A.leading_dimension(), A.buffer().get_stream().
                                            raw_stream()));
   A.buffer().synchronization_point();
}

inline
void
gemm(cublas::handle& H, char Atrans, char Btrans, int M, int N, int K, double alpha,
     cuda::gpu_buffer<double> const& A, int lda, cuda::gpu_buffer<double> const& B, int ldb,
     double beta, cuda::gpu_buffer<double>& C, int ldc)
{
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   check_error(cublasDgemm(H.raw_handle(), cublasOperation_t(Atrans), cublasOperation_t(Btrans), M, N, K,
                           &alpha, A.device_ptr(), lda, B.device_ptr(), ldb,
                           &beta, C.device_ptr(), ldc));
   C.synchronization_point();
}

template <typename T, typename U, typename V>
inline
void
gemm(T alpha, blas::BlasMatrix<T, gpu_matrix<T>, U> const& A, blas::BlasMatrix<T, gpu_matrix<T>, V> const& B,
     T beta, gpu_matrix<T>& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   gemm(get_handle(), A.trans(), B.trans(), A.rows(), A.cols(), B.cols(), alpha, A.as_derived().buffer(),
	A.leading_dimension(), B.as_derived().buffer(), B.leading_dimension(), beta,
        C.buffer(), C.leading_dimension());
}

} // namespace cublas

#endif
