// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cublas/cub.h
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

// Miscellaneous GPU operations using the CUB library

#if !defined(MPTOOLKIT_CUDA_CUB_H)
#define MPTOOLKIT_CUDA_CUB_H

#include "cuda.h"
#include "gpu_buffer.h"

namespace cuda
{

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx, cuda::gpu_ref<T>& r);

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx, cuda::gpu_ref<T>&& r);

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx, cuda::gpu_ref<std::complex<T>>& r);

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx, cuda::gpu_ref<std::complex<T>>&& r);

template <typename T>
void vector_fill(T alpha, int N, cuda::gpu_ptr<T> y, int incy);

template <typename T>
void
vector_permute(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::gpu_ptr<T> y, int incy, int const* Perm);

template <typename T>
void
vector_parallel(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::const_gpu_ptr<T> y, int incy,
		cuda::gpu_ptr<T> z, int incz);

template <typename T>
void
vector_add_parallel(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::const_gpu_ptr<T> y, int incy,
		    cuda::gpu_ptr<T> z, int incz);

template <typename T>
void
vector_parallel_scaled(int n, T const& alpha, cuda::const_gpu_ptr<T> x, int incx, cuda::const_gpu_ptr<T> y, int incy,
		       cuda::gpu_ptr<T> z, int incz);

template <typename T>
void
vector_add_parallel_scaled(int n, T const& alpha, cuda::const_gpu_ptr<T> x, int incx, cuda::const_gpu_ptr<T> y, int incy,
			   cuda::gpu_ptr<T> z, int incz);

template <typename T>
void matrix_fill(T alpha, int M, int N, cuda::gpu_ptr<T> A, int lda);

void
dgmm(int M, int K,
     cuda::const_gpu_ptr<double> x, int incx,
     cuda::const_gpu_ptr<std::complex<double>> B, int ldb,
     cuda::gpu_ptr<std::complex<double>> C, int ldc);


} // namespace cuda

#endif
