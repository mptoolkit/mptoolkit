// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cublas/gpu_vector.h
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

#if !defined(MPTOOLKIT_CUBLAS_GPU_VECTOR_H)
#define MPTOOLKIT_CUBLAS_GPU_VECTOR_H

#include "cublas.h"
#include "gpu_buffer.h"
#include "blas/matrixref.h"
#include "blas/matrix.h"

namespace cublas
{

//
// GPU-storage vector type.
//
template <typename T>
class gpu_vector;

} // namespace cublas

namespace blas
{
template <typename T>
struct blas_traits<cublas::gpu_vector<T>>
{
   using storage_type       = cuda::gpu_ptr<T>;
   using const_storage_type = cuda::const_gpu_ptr<T>;
};

} // namespace blas

namespace cublas
{

template <typename T>
class gpu_vector : public blas::BlasVector<T, gpu_vector<T>>
{
   public:
      using value_type         = T;
      using storage_type       = cuda::gpu_ptr<T>;
      using const_storage_type = cuda::const_gpu_ptr<T>;

      gpu_vector(int Size_);

      gpu_vector(int Size_, blas::arena const& A);

      gpu_vector(gpu_vector&& Other) = default;

      gpu_vector(gpu_vector const&) = delete;

      gpu_vector& operator=(gpu_vector&&) = delete;

      ~gpu_vector()
      {
      }

      gpu_vector& operator=(gpu_vector const& Other)
      {
         assign(*this, Other);
         return *this;
      }

      // assignment of expressions based on the same vector type
      template <typename U>
      gpu_vector& operator=(blas::VectorRef<T, gpu_vector<T>, U> const& E)
      {
	 assign(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      gpu_vector& operator+=(blas::VectorRef<T, gpu_vector<T>, U> const& E)
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      gpu_vector& operator-=(blas::VectorRef<T, gpu_vector<T>, U> const& E)
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      constexpr int stride() const { return 1; }

      int size() const { return Size; }

      cuda::gpu_buffer<T>& buffer() { return Buf; }
      cuda::gpu_buffer<T> const& buffer() const { return Buf; }

      storage_type storage() { return Buf.ptr(); }
      const_storage_type storage() const { return Buf.cptr(); }

   private:
      int Size;
      cuda::gpu_buffer<T> Buf;
};

} // namespace cublas

namespace blas
{

template <typename T>
struct default_arena<cublas::gpu_vector<T>> : cublas::detail::gpu_default_arena
{
   // Arena is inhereted from the base class
};

} // namespace blas

namespace cublas
{

template <typename T>
inline
gpu_vector<T>::gpu_vector(int Size_, blas::arena const& Arena)
   : Size(Size_), Buf(cuda::gpu_buffer<T>::allocate(Size_, Arena))
{
}

template <typename T>
inline
gpu_vector<T>::gpu_vector(int Size_)
   : gpu_vector(Size_, blas::default_arena<gpu_vector<T>>::Arena)
{
}

// blocking vector get
template <typename T>
blas::Vector<T>
get_wait(gpu_vector<T> const& M)
{
   blas::Vector<T> Result(M.size());
   cublas::check_error(cublasGetVector(M.size(), sizeof(T),
				       M.storage().device_ptr(), M.stride(),
				       Result.storage(), Result.stride()));
   return Result;
}

// blocking vector set
template <typename T>
void
set_wait(gpu_vector<T>& A, blas::Vector<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   //TRACE(A.cols())(A.rows())(A.device_ptr())(A.leading_dim())(B.data())(leading_dimension(B));
   cublas::check_error(cublasSetVector(A.size(), sizeof(T),
				       B.storage().device_ptr(), B.stride(),
				       A.storage().device_ptr(), A.stride()));
}

// non-blocking set
template <typename T>
cuda::event
set(gpu_vector<T>& A, blas::Vector<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.size(), B.size());
   cublas::setVectorAsync(A.size(), B.storage(), B.stride(),
                          A.storage(), A.stride());
   return A.storage().sync();
}

// copy

template <typename T>
gpu_vector<T>
copy(gpu_vector<T> const& x, blas::arena const& A)
{
   gpu_vector<T> Result(x.size(), A);
   Result = x;
   return Result;
}

// BLAS functions

template <typename T, typename U, typename V>
inline
void
gemv(T alpha, blas::BlasMatrix<T, gpu_matrix<T>, U> const& A,
     blas::BlasVector<T, gpu_vector<T>, V> const& x, T beta,
     gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   cublas::gemv(get_handle(), A.trans(), A.rows(), A.cols(), alpha, A.storage(),
                A.leading_dimension(), x.storage(), x.stride(), beta,
                y.storage(), y.stride());
}

template <typename T, typename U>
inline
void vector_copy_scaled(T alpha, blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   cublas::copy(get_handle(), x.size(), x.storage(), x.stride(), y.storage(), y.stride());
   cublas::scal(get_handle(), y.size(), alpha, y.storage(), y.stride());
}

template <typename T, typename U>
inline
void vector_copy(blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   cublas::copy(get_handle(), x.size(), x.storage(), x.stride(), y.storage(), y.stride());
}

template <typename T, typename U>
inline
void vector_add_scaled(T alpha, blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   cublas::axpy(get_handle(), x.size(), alpha,
                x.storage(), x.stride(),
                y.storage(), y.stride());
}

template <typename T, typename U>
inline
void vector_add(blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   cublas::axpy(get_handle(), x.size(), blas::number_traits<T>::identity(),
                x.storage(), x.stride(),
                y.storage(), y.stride());
}

template <typename T>
inline
void vector_scale(T alpha, gpu_vector<T>& y)
{
   cublas::scal(get_handle(), y.size(), alpha, y.storage(), y.stride());
}


} // namespace cublas

#endif
