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

#include "blas/matrixref.h"
#include "blas/matrix.h"
#include "blas/vector_view.h"
#include "cublas.h"
#include "gpu_vector.h"
#include "gpu_buffer.h"

namespace cublas
{

//
// GPU-storage matrix type.  Column-major format for compatability with BLAS
//

template <typename T>
using gpu_vector_view = blas::vector_view<T, gpu_tag>;

template <typename T>
using const_gpu_vector_view = blas::const_vector_view<T, gpu_tag>;

template <typename T>
class gpu_matrix : public blas::NormalMatrix<T, gpu_matrix<T>, gpu_tag>
{
   public:
      using value_type         = T;
      using storage_type       = cuda::gpu_ptr<T>;
      using const_storage_type = cuda::const_gpu_ptr<T>;

      gpu_matrix(int Rows_, int Cols_);

      gpu_matrix(int Rows_, int Cols_, int leadingdim);

      gpu_matrix(int Rows_, int Cols_, blas::arena const& A, int leadingdim);

      gpu_matrix(int Rows_, int Cols_, blas::arena const& A);

      gpu_matrix(gpu_matrix&& Other) = default;

      gpu_matrix(gpu_matrix const&) = delete;

      template <typename U>
      gpu_matrix(blas::MatrixRef<T, U, gpu_tag> const& E)
         : gpu_matrix(E.rows(), E.cols())
      {
         assign(*this, E.as_derived());
      }

      gpu_matrix& operator=(gpu_matrix&&) = delete;

      ~gpu_matrix() = default;

      gpu_matrix& operator=(gpu_matrix const& Other)
      {
         assign(*this, Other);
         return *this;
      }

      // assignment of expressions based on the same matrix type
      template <typename U>
      gpu_matrix& operator=(blas::MatrixRef<T, U, gpu_tag> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      gpu_matrix& operator+=(blas::MatrixRef<T, U, gpu_tag> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      gpu_matrix& operator-=(blas::MatrixRef<T, U, gpu_tag> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      constexpr char trans() const { return 'N'; }

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dimension() const { return LeadingDimension; }

      gpu_vector_view<T>
      row(int r)
      {
         return gpu_vector_view<T>(Cols, LeadingDimension, Buf.ptr(r));
      }

      const_gpu_vector_view<T>
      row(int r) const
      {
         return const_gpu_vector_view<T>(Cols, LeadingDimension, Buf.ptr(r));
      }

      gpu_vector_view<T>
      column(int c)
      {
         return gpu_vector_view<T>(Rows, 1, Buf.ptr(LeadingDimension*c));
      }

      const_gpu_vector_view<T>
      column(int c) const
      {
         return const_gpu_vector_view<T>(Rows, 1, Buf.ptr(LeadingDimension*c));
      }

      gpu_vector_view<T>
      diagonal()
      {
         return gpu_vector_view<T>(std::min(Rows,Cols), LeadingDimension+1, Buf.ptr());
      }

      const_gpu_vector_view<T>
      diagonal() const
      {
         return const_gpu_vector_view<T>(std::min(Rows,Cols), LeadingDimension+1, Buf.ptr());
      }

      cuda::gpu_buffer<T>& buffer() { return Buf; }
      cuda::gpu_buffer<T> const& buffer() const { return Buf; }

      storage_type storage() { return Buf.ptr(); }
      const_storage_type storage() const { return Buf.cptr(); }

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

} // namespace cublas

namespace blas
{

template <typename T>
struct default_arena<cublas::gpu_matrix<T>> : cublas::detail::gpu_default_arena
{
   // Arena is inhereted from the base class
};

} // namespace blas

namespace cublas
{

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, blas::arena const& Arena, int leadingdim)
   : Rows(Rows_), Cols(Cols_), LeadingDimension(leadingdim),
     Buf(cuda::gpu_buffer<T>::allocate(LeadingDimension*Rows, Arena))
{
   DEBUG_CHECK(LeadingDimension >= Rows);
}

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, blas::arena const& Arena)
   : gpu_matrix(Rows_, Cols_, Arena, select_leading_dimension(Rows_))
{
}

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, int leadingdim)
   : gpu_matrix(Rows_, Cols_, blas::default_arena<gpu_matrix<T>>::Arena, leadingdim)
{
}

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_)
   : gpu_matrix(Rows_, Cols_, blas::default_arena<gpu_matrix<T>>::Arena, select_leading_dimension(Rows_))
{
}

// blocking matrix get
template <typename T>
blas::Matrix<T>
get_wait(gpu_matrix<T> const& M)
{
   blas::Matrix<T> Result(M.rows(), M.cols());
   cublas::check_error(cublasGetMatrix(M.rows(), M.cols(), sizeof(T),
				       M.storage().device_ptr(), M.leading_dimension(),
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
                          A.storage(), A.leading_dimension());
   return A.storage().sync();
}

// copy

template <typename T>
inline
gpu_matrix<T>
copy(gpu_matrix<T> const& x, blas::arena const& A)
{
   gpu_matrix<T> Result(x.rows(), x.cols(), A, x.leading_dimension());
   Result = x;
   return Result;
}

} // namespace cublas

#endif
