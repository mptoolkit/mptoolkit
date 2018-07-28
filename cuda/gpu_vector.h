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
#include "blas/vector_view.h"

namespace blas
{

//
// GPU-storage vector type.
//

// defined in cublas.cpp
namespace detail
{
extern arena gpu_default_arena;
} // namespace detail

struct gpu_tag
{
   template <typename T>
   using buffer_type = cuda::gpu_buffer<T>;

   template <typename T>
   using storage_type = cuda::gpu_ptr<T>;

   template <typename T>
   using const_storage_type = cuda::const_gpu_ptr<T>;

   template <typename T>
   using async_proxy = cuda::gpu_ref<T>;

   template <typename T>
   using async_ref = cuda::gpu_ref<T>;

   template <typename T>
   static
   inline
   async_ref<T> allocate_async_ref()
   {
      return cuda::allocate_gpu_ref<T>();
   }

   template <typename T>
   static
   blas::arena default_arena() { return detail::gpu_default_arena; }

   template <typename T>
   static int select_leading_dimension(int ld)
   {
      return  ld == 1 ? 1 : cuda::round_up(ld, 32);
   }

   template <typename T>
   static void uninitialized_default_construct_n(cuda::gpu_ptr<T> p, int n);

   template <typename T>
   static void uninitialized_fill_n(cuda::gpu_ptr<T> p, int n, T fill);

   template <typename T>
   static void destroy_n(cuda::gpu_ptr<T> p, int n);
};

template <typename T>
inline
void gpu_tag::uninitialized_default_construct_n(cuda::gpu_ptr<T> p, int n)
{
   // cuda types have trivial constructors (or std::complex, we assume
   // doesn't need the constructor to be called)
   static_assert(cuda::is_cuda_floating_point_v<T>);
}

template <typename T>
inline
void gpu_tag::uninitialized_fill_n(cuda::gpu_ptr<T> p, int n, T fill)
{
   // cuda types have trivial constructors (or std::complex, we assume
   // doesn't need the constructor to be called)
   static_assert(cuda::is_cuda_floating_point_v<T>);
   vector_fill(fill, n, p, 1);
}

template <typename T>
inline
void gpu_tag::destroy_n(cuda::gpu_ptr<T> p, int n)
{
   // cuda types have trivial constructors (or std::complex, we assume
   // doesn't need the constructor to be called)
   static_assert(cuda::is_cuda_floating_point_v<T>);
}

template <typename T>
using gpu_vector = Vector<T, gpu_tag>;

// blocking vector get
template <typename T, typename U>
blas::Vector<T>
get_wait(blas::BlasVector<T, U, gpu_tag> const& M)
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
set_wait(gpu_vector<T>& A, blas::BlasVector<T, U, blas::cpu_tag> const& B)
{
   DEBUG_CHECK_EQUAL(A.size(), B.size());
   cublas::setVector(A.size(), B.storage(), B.stride(), A.storage(), A.stride());
}

template <typename T, typename U, typename V>
void
set_wait(blas::BlasVectorProxy<T, U, gpu_tag>&& A, blas::BlasVector<T, V, blas::cpu_tag> const& B)
{
   DEBUG_CHECK_EQUAL(A.size(), B.size());
   cublas::setVector(A.size(), B.storage(), B.stride(), std::move(A).storage(), A.stride());
}

// non-blocking set
template <typename T, typename U>
cuda::event
set(gpu_vector<T>& A, blas::BlasVector<T, U, blas::cpu_tag> const& B)
{
   DEBUG_CHECK_EQUAL(A.size(), B.size());
   cublas::setVectorAsync(A.size(), B.storage(), B.stride(),
                          A.storage(), A.stride());
   return A.storage().sync();
}

template <typename T, typename U, typename V>
cuda::event
set(blas::BlasVectorProxy<T, U, gpu_tag>&& A, blas::BlasVector<T, V, blas::cpu_tag> const& B)
{
   DEBUG_CHECK_EQUAL(A.size(), B.size());
   cublas::setVectorAsync(A.size(), B.storage(), B.stride(),
                          std::move(A).storage(), A.stride());
   return A.storage().sync();
}

} // namespace blas

#endif
