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
class gpu_matrix;

} // namespace cublas

namespace blas
{
template <typename T>
struct blas_traits<cublas::gpu_matrix<T>>
{
   using storage_type       = cuda::gpu_buffer<T>*;
   using const_storage_type = cuda::gpu_buffer<T> const*;
};

} // namespace blas

namespace cublas
{

template <typename T>
class gpu_matrix : public blas::BlasMatrix<T, gpu_matrix<T>>
{
   public:
      using value_type         = T;
      using storage_type       = cuda::gpu_buffer<T>*;
      using const_storage_type = cuda::gpu_buffer<T> const*;

      gpu_matrix(int Rows_, int Cols_, arena const& A, int leadingdim);

      gpu_matrix(int Rows_, int Cols_, arena const& A);

      gpu_matrix(gpu_matrix&& Other) = default;

      gpu_matrix(gpu_matrix const&) = delete;

      gpu_matrix& operator=(gpu_matrix&&) = delete;

      ~gpu_matrix()
      {
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

      template <typename U>
      gpu_matrix& operator+=(blas::MatrixRef<T, gpu_matrix<T>, U> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      gpu_matrix& operator-=(blas::MatrixRef<T, gpu_matrix<T>, U> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      constexpr char trans() const { return 'N'; }

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dimension() const { return LeadingDimension; }

      cuda::gpu_buffer<T>& buffer() { return Buf; }
      cuda::gpu_buffer<T> const& buffer() const { return Buf; }

      storage_type storage() { return &Buf; }
      const_storage_type storage() const { return &Buf; }

      static int select_leading_dimension(int ld)
      {
         return ld == 1 ? 1 : cuda::round_up(ld, 32);
      }

   private:
      int Rows;
      int Cols;
      int LeadingDimension;
      cuda::gpu_buffer<T> Buf;
};

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, arena const& Arena, int leadingdim)
   : Rows(Rows_), Cols(Cols_), LeadingDimension(leadingdim),
     Buf(cuda::gpu_buffer<T>::allocate(LeadingDimension*Rows, Arena))
{
   DEBUG_CHECK(LeadingDimension >= Rows);
}

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, arena const& Arena)
   : gpu_matrix(Rows_, Cols_, Arena, select_leading_dimension(Rows_))
{
}

// blocking matrix get
template <typename T>
blas::Matrix<T>
get_wait(gpu_matrix<T> const& M)
{
   blas::Matrix<T> Result(M.rows(), M.cols());
   cublas::check_error(cublasGetMatrix(M.rows(), M.cols(), sizeof(T),
				       M.storage()->device_ptr(), M.leading_dimension(),
				       Result.storage(), Result.leading_dimension()));
   return Result;
}

// blocking matrix set
template <typename T>
void
set_wait(gpu_matrix<T>& A, blas::Matrix<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   //TRACE(A.cols())(A.rows())(A.device_ptr())(A.leading_dim())(B.data())(leading_dimension(B));
   cublas::check_error(cublasSetMatrix(A.rows(), A.cols(), sizeof(T),
				       B.data(), leading_dimension(B),
				       A.device_ptr(), A.leading_dim()));
}

// non-blocking set
template <typename T>
cuda::event
set(gpu_matrix<T>& A, blas::Matrix<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   cublas::setMatrixAsync(A.rows(), A.cols(), B.storage(), B.leading_dimension(),
                          *A.storage(), A.leading_dimension());
   return A.storage()->sync();
}

// BLAS functions

template <typename T, typename U, typename V>
inline
void
gemm(T alpha, blas::BlasMatrix<T, gpu_matrix<T>, U> const& A,
     T beta, blas::BlasMatrix<T, gpu_matrix<T>, V> const& B,
     gpu_matrix<T>& C)
{
   DEBUG_CHECK_EQUAL(A.cols(), B.rows());
   DEBUG_CHECK_EQUAL(A.rows(), C.rows());
   DEBUG_CHECK_EQUAL(B.cols(), C.cols());
   cublas::gemm(get_handle(), A.trans(), B.trans(), A.rows(), A.cols(), B.cols(), alpha, *A.storage(),
                A.leading_dimension(), *B.storage(), B.leading_dimension(), beta,
                *C.storage(), C.leading_dimension());
}

template <typename T, typename U>
inline
void matrix_copy_scaled(T alpha, blas::BlasMatrix<T, gpu_matrix<T>, U> const& A, gpu_matrix<T>& C)
{
   cublas::geam(get_handle(), A.trans(), A.rows(), A.cols(),
                alpha, A.storage(), A.leading_dimension(),
                blas::number_traits<T>::zero(), C.storage(), C.leading_dimension());
}

template <typename T, typename U>
inline
void matrix_copy(blas::BlasMatrix<T, gpu_matrix<T>, U> const& A, gpu_matrix<T>& C)
{
   cublas::geam(get_handle(), A.trans(), A.rows(), A.cols(),
                blas::number_traits<T>::identity(), A.storage(), A.leading_dimension(),
                blas::number_traits<T>::zero(), C.storage(), C.leading_dimension());
}

template <typename T, typename U>
inline
void matrix_add_scaled(T alpha, blas::BlasMatrix<T, gpu_matrix<T>, U> const& A, gpu_matrix<T>& C)
{
   cublas::geam(get_handle(), A.trans(), A.rows(), A.cols(),
                alpha, A.storage(), A.leading_dimension(),
                blas::number_traits<T>::identity(), C.storage(), C.leading_dimension());
}

template <typename T, typename U>
inline
void matrix_add(blas::BlasMatrix<T, gpu_matrix<T>, U> const& A, gpu_matrix<T>& C)
{
   cublas::geam(get_handle(), A.trans(), A.rows(), A.cols(),
                blas::number_traits<T>::identity(), A.storage(), A.leading_dimension(),
                blas::number_traits<T>::identity(), C.storage(), C.leading_dimension());
}

template <typename T, typename U, typename V>
inline
void geam(T alpha, blas::BlasMatrix<T, gpu_matrix<T>, U> const& A,
          T beta, blas::BlasMatrix<T, gpu_matrix<T>, U> const& B,
          gpu_matrix<T>& C)
{
   cublas::geam(get_handle(), A.trans(), B.trans(), A.rows(), A.cols(),
                alpha, A.storage(), A.leading_dimension(),
                beta, B.storage(), B.leading_dimension(),
                C.storage(), C.leading_dimension());
}

} // namespace cublas

#endif
