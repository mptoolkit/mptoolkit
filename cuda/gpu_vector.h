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
#include "blas/vectorref.h"
#include "blas/vector.h"

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

      gpu_vector() = delete;

      gpu_vector(int Size_);

      gpu_vector(int Size_, blas::arena const& A);

      gpu_vector(gpu_vector&& Other) = default;

      gpu_vector(gpu_vector const&) = delete;

      gpu_vector& operator=(gpu_vector&&) = delete;

      ~gpu_vector() = default;

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

//
// gpu_vector_view is a proxy class that interprets strided 'view' of a
// gpu memory buffer as a vector.  It can be used as an l-value.
//

template <typename T>
class gpu_vector_view : public blas::BlasVector<T, gpu_vector<T>, gpu_vector_view<T>>
{
   public:
      using value_type         = T;
      using storage_type       = cuda::gpu_ptr<T>;
      using const_storage_type = cuda::const_gpu_ptr<T>;

      gpu_vector_view() = delete;

      gpu_vector_view(int Size_, int Stride_, cuda::gpu_ptr<T> Ptr_)
         : Size(Size_), Stride(Stride_), Ptr(Ptr_) {}

      gpu_vector_view(gpu_vector_view&& Other) = default;

      gpu_vector_view(gpu_vector_view const&) = delete;

      gpu_vector_view& operator=(gpu_vector_view&&) = delete;

      ~gpu_vector_view() = default;

      template <typename U>
      gpu_vector_view&& operator=(blas::VectorRef<T, gpu_vector<T>, U> const& E) &&
      {
	 assign(std::move(*this), E.as_derived());
	 return std::move(*this);
      }

      template <typename U>
      gpu_vector_view&& operator+=(blas::VectorRef<T, gpu_vector<T>, U> const& E) &&
      {
	 add(*this, E.as_derived());
	 return *this;
      }

      template <typename U>
      gpu_vector_view&& operator-=(blas::VectorRef<T, gpu_vector<T>, U> const& E) &&
      {
	 subtract(*this, E.as_derived());
	 return *this;
      }

      int stride() const { return Stride; }

      int size() const { return Size; }

      storage_type storage() && { return Ptr; }
      const_storage_type storage() const& { return Ptr; }

   private:
      int Size;
      int Stride;
      cuda::gpu_ptr<T> Ptr;
};

template <typename T>
class const_gpu_vector_view : public blas::BlasVector<T, gpu_vector<T>, gpu_vector_view<T>>
{
   public:
      using value_type         = T;
      using storage_type       = cuda::const_gpu_ptr<T>;
      using const_storage_type = cuda::const_gpu_ptr<T>;

      const_gpu_vector_view() = delete;

      const_gpu_vector_view(int Size_, int Stride_, cuda::const_gpu_ptr<T> Ptr_);

      const_gpu_vector_view(const_gpu_vector_view&& Other) = default;

      const_gpu_vector_view(const_gpu_vector_view const&) = delete;

      const_gpu_vector_view& operator=(const_gpu_vector_view&&) = delete;

      ~const_gpu_vector_view() = default;

      int stride() const { return Stride; }

      int size() const { return Size; }

      const_storage_type storage() const { return Ptr; }

   private:
      int Size;
      int Stride;
      cuda::const_gpu_ptr<T> Ptr;
};

namespace detail
{
struct gpu_default_arena
{
   static blas::arena Arena;
};
} // namespace detail

namespace detail
{
blas::arena gpu_default_arena::Arena = blas::arena(new cuda::BlockAllocator(cuda::DefaultBlockMultiple, false));
} // namespace detail

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
template <typename T, typename U>
blas::Vector<T>
get_wait(blas::BlasVector<T, gpu_vector<T>, U> const& M)
{
   blas::Vector<T> Result(M.size());
   cublas::check_error(cublasGetVector(M.size(), sizeof(T),
				       M.storage().device_ptr(), M.stride(),
				       Result.storage(), Result.stride()));
   return Result;
}

// blocking vector set
template <typename T, typename U>
void
set_wait(gpu_vector<T>& A, blas::BlasVector<T, blas::Vector<T>, U> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   cublas::setVector(A.size(), B.storage(), B.stride(), A.storage(), A.stride());
}

template <typename T, typename U>
void
set_wait(gpu_vector_view<T>&& A, blas::BlasVector<T, blas::Vector<T>, U> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   cublas::setVector(A.size(), B.storage(), B.stride(), A.storage(), A.stride());
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

template <typename T>
cuda::event
set(gpu_vector_view<T>&& A, blas::Vector<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.size(), B.size());
   cublas::setVectorAsync(A.size(), B.storage(), B.stride(),
                          A.storage(), A.stride());
   return A.storage().sync();
}

// copy

template <typename T, typename U>
gpu_vector<T>
copy(blas::BlasVector<T, gpu_vector<T>, U> const& x, blas::arena const& A)
{
   gpu_vector<T> Result(x.size(), A);
   assign(Result, x.derived());
   return Result;
}

template <typename T, typename U>
gpu_vector<T>
copy(blas::BlasVector<T, gpu_vector<T>, U> const& x)
{
   gpu_vector<T> Result(x.size());
   assign(Result, x.derived());
   return Result;
}

// BLAS functions

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
void vector_copy_scaled(T alpha, blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector_view<T>&& y)
{
   DEBUG_CHECK_EQUAL(x.size(), y.size());
   cublas::copy(get_handle(), x.size(), x.storage(), x.stride(), std::move(y).storage(), y.stride());
   cublas::scal(get_handle(), y.size(), alpha, std::move(y).storage(), y.stride());
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
void vector_copy(blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector_view<T>&& y)
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
void vector_add_scaled(T alpha, blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector_view<T>&& y)
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

template <typename T, typename U>
inline
void vector_add(blas::BlasVector<T, gpu_vector<T>, U> const& x, gpu_vector_view<T>&& y)
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

template <typename T>
inline
void vector_scale(T alpha, gpu_vector_view<T>&& y)
{
   cublas::scal(get_handle(), y.size(), alpha, y.storage(), y.stride());
}

} // namespace cublas

#endif
