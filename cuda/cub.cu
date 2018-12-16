// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cublas/cub.cu
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

#include "cub/cub.cuh"
#include "cuda/cub.h"
#include "cuda/gpu_buffer.h"
#include <type_traits>
#include "common/trace.h"

namespace cuda
{

template <typename T>
__global__ void cuda_vector_fill(T alpha, unsigned n, T* y, unsigned incy)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      y[i*incy] = alpha;
}

template <typename T>
void vector_fill(T alpha, int n, cuda::gpu_ptr<T> y, int incy)
{
   unsigned threads_per_block = 64;
   unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
   cuda_vector_fill<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>(alpha, n, y.device_ptr(), incy);
}

// template instantiations
template void vector_fill<float>(float alpha, int n, cuda::gpu_ptr<float> y, int incy);
template void vector_fill<double>(double alpha, int n, cuda::gpu_ptr<double> y, int incy);
template void vector_fill<std::complex<float>>(std::complex<float> alpha, int n,
					       cuda::gpu_ptr<std::complex<float>> y, int incy);
template void vector_fill<std::complex<double>>(std::complex<double> alpha, int n,
						cuda::gpu_ptr<std::complex<double>> y, int incy);


template <typename T>
__global__ void cuda_matrix_fill(T alpha, unsigned M, unsigned N, T* y, unsigned lda)
{
   unsigned Col = blockIdx.x*blockDim.x + threadIdx.x;
   unsigned Row = blockIdx.y*blockDim.y + threadIdx.y;
   if (Row < M && Col < N)
      y[lda*Col + Row] = alpha;
}

template <typename T>
void matrix_fill(T alpha, int M, int N, cuda::gpu_ptr<T> A, int lda)
{
   dim3 BlockDim(16, 16);
   dim3 GridDim((M+BlockDim.x-1) / BlockDim.x, (N+BlockDim.y-1) / BlockDim.y);
   cuda_matrix_fill<<<GridDim, BlockDim, 0, A.get_stream().raw_stream()>>>(alpha, M, N, A.device_ptr(), lda);
}

// template instantiations
template void matrix_fill<float>(float alpha, int m, int n, cuda::gpu_ptr<float> y, int lda);
template void matrix_fill<double>(double alpha, int m, int n, cuda::gpu_ptr<double> y, int lda);
template void matrix_fill<std::complex<float>>(std::complex<float> alpha, int m, int n,
					       cuda::gpu_ptr<std::complex<float>> y, int lda);
template void matrix_fill<std::complex<double>>(std::complex<double> alpha, int m, int n,
						cuda::gpu_ptr<std::complex<double>> y, int lda);



namespace detail
{

template <typename T>
class stride_ptr
{
   public:
      using value_type        = typename std::remove_const<T>::type;
      using iterator_category = std::random_access_iterator_tag;
      using pointer           = T*;
      using reference         = T&;
      using difference_type   = std::ptrdiff_t;

      // CUDA 9 compiler doesn't like __host__ or __device__ designations on trivial constructors
      // __host__ __device__ stride_ptr() {}
      // __host__ __device__ stride_ptr(stride_ptr const& Other) = default;

      __host__ __device__ stride_ptr(T* x, int stride) : x_(x), stride_(stride) {}

      __host__ __device__ reference operator*() const { return *x_; }
      __host__ __device__ pointer operator->() const { return x_; }

      __host__ __device__ reference operator[](int i) const { return x_[i*stride_]; }

      __host__ __device__ stride_ptr& operator++() { x_ += stride_; return this; }
      __host__ __device__ stride_ptr operator++(int) { T* Temp = x_; x_ += stride_; return stride_ptr(Temp, stride_); }

      __host__ __device__ stride_ptr& operator--() { x_ -= stride_; return this; }
      __host__ __device__ stride_ptr operator--(int) { T* Temp = x_; x_ -= stride_; return stride_ptr(Temp, stride_); }

      __host__ __device__ stride_ptr& operator+=(int i) { x_ += i*stride_; return *this; }

      __host__ __device__ int stride() const { return stride_; }

   private:
      T* x_;
      int stride_;
};

template <typename T>
inline
__host__ __device__
stride_ptr<T>
operator+(stride_ptr<T> const& x, int i)
{
   return stride_ptr<T>(x.operator->()+i*x.stride(), x.stride());
}

} // namespace detail

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx, cuda::gpu_ref<T>& r)
{
   TRACE_CUDA("vector_sum")(Size)(x)(incx)(r);
   std::size_t TempStorageBytes = 0;
   cuda::check_error(cub::DeviceReduce::Sum(nullptr, TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
					    r.device_ptr(), Size));
   if (TempStorageBytes == 0)
      TempStorageBytes = 1;
   gpu_buffer<unsigned char> TempBuffer = allocate_gpu_temporary<unsigned char>(TempStorageBytes);
   TempBuffer.wait_for(x);
   TempBuffer.wait_for(r);
   cuda::check_error(cub::DeviceReduce::Sum(static_cast<void*>(TempBuffer.device_ptr()),
					    TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
					    r.device_ptr(), Size, TempBuffer.get_stream().raw_stream()));
   r.wait(TempBuffer.sync());
   x.wait(r.sync());
}

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx,
	   cuda::gpu_ref<std::complex<T>>& r)
{
   TRACE_CUDA("vector_sum")(Size)(x)(incx)(r);
   std::size_t TempStorageBytes1 = 0;
   // real part
   // workspace query
   cuda::check_error(cub::DeviceReduce::Sum(nullptr, TempStorageBytes1,
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
					    static_cast<T*>(static_cast<void*>(r.device_ptr())), Size));
   if (TempStorageBytes1 == 0)
      TempStorageBytes1 = 1;
   gpu_buffer<unsigned char> TempBuffer1 = allocate_gpu_temporary<unsigned char>(TempStorageBytes1);
   TempBuffer1.wait_for(x);
   TempBuffer1.wait_for(r);
   // do the work
   cuda::check_error(cub::DeviceReduce::Sum(static_cast<void*>(TempBuffer1.device_ptr()), TempStorageBytes1,
					    detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
					    static_cast<T*>(static_cast<void*>(r.device_ptr())), Size,
					    TempBuffer1.get_stream().raw_stream()));
   //   cuda::device_synchronize();

   // imag part
   // workspace query
   std::size_t TempStorageBytes2 = 0;
   cuda::check_error(cub::DeviceReduce::Sum(nullptr, TempStorageBytes2,
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
					    static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size));
   if (TempStorageBytes2 == 0)
      TempStorageBytes2 = 1;
   gpu_buffer<unsigned char> TempBuffer2 = allocate_gpu_temporary<unsigned char>(TempStorageBytes2);
   TempBuffer2.wait_for(x);
   TempBuffer2.wait_for(r);
   // do the work
   cuda::check_error(cub::DeviceReduce::Sum(static_cast<void*>(TempBuffer2.device_ptr()), TempStorageBytes2,
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size,
					    TempBuffer2.get_stream().raw_stream()));

   //   cuda::device_synchronize();

   r.wait(TempBuffer1.sync());
   r.wait(TempBuffer2.sync());
   x.wait(r.sync());
}

// template instantiations
template void vector_sum<double>(int Size, cuda::const_gpu_ptr<double> const& x, int incx,
				 cuda::gpu_ref<double>& r);

template void vector_sum<double>(int Size,
				 cuda::const_gpu_ptr<std::complex<double>> const& x, int incx,
				 cuda::gpu_ref<std::complex<double>>& r);

template void vector_sum<float>(int Size, cuda::const_gpu_ptr<float> const& x, int incx,
				cuda::gpu_ref<float>& r);

template void vector_sum<float>(int Size,
				cuda::const_gpu_ptr<std::complex<float>> const& x, int incx,
				cuda::gpu_ref<std::complex<float>>& r);

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx,
	   cuda::gpu_ref<std::complex<T>>&& r)
{
   vector_sum(Size, x, incx, r);
}

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx,
	   cuda::gpu_ref<T>&& r)
{
   vector_sum(Size, x, incx, r);
}

template void vector_sum<double>(int Size, cuda::const_gpu_ptr<double> const& x, int incx,
				 cuda::gpu_ref<double>&& r);

template void vector_sum<double>(int Size,
				 cuda::const_gpu_ptr<std::complex<double>> const& x, int incx,
				 cuda::gpu_ref<std::complex<double>>&& r);

template void vector_sum<float>(int Size, cuda::const_gpu_ptr<float> const& x, int incx,
				cuda::gpu_ref<float>&& r);

template void vector_sum<float>(int Size,
				cuda::const_gpu_ptr<std::complex<float>> const& x, int incx,
				cuda::gpu_ref<std::complex<float>>&& r);

template <typename T>
__global__ void permute_vector_cuda(unsigned n, T const* In, T* Out, int const* Perm)
{
   unsigned i =  blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      Out[i] = In[Perm[i]];
}

template <typename T>
__global__ void permute_vector_stride_cuda(unsigned n, T const* In, int InStride, T* Out, int OutStride,
						  int const* Perm)
{
   unsigned i =  blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      Out[i*OutStride] = In[Perm[i]*InStride];
}

template <typename T>
void
vector_permute(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::gpu_ptr<T> y, int incy, int const* Perm)
{
   y.wait_for(x);
   gpu_buffer<int> PermDevice = allocate_gpu_temporary<int>(n);
   memcpy_host_to_device_async(PermDevice.get_stream(), Perm, PermDevice.device_ptr(), n*sizeof(int));
   PermDevice.wait_for(y);
   PermDevice.wait_for(x);
   if (incx == 1 && incy == 1)
   {
      unsigned threads_per_block = 64;
      unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
      permute_vector_cuda<<<nblocks,threads_per_block, 0,
	 PermDevice.get_stream().raw_stream()>>>(n, x.device_ptr(), y.device_ptr(),
						 PermDevice.device_ptr());
   }
   else
   {
      unsigned threads_per_block = 64;
      unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
      permute_vector_stride_cuda<<<nblocks,threads_per_block, 0,
	 PermDevice.get_stream().raw_stream()>>>(n, x.device_ptr(), incx,
						 y.device_ptr(), incy,
						 PermDevice.device_ptr());
   }
   x.wait_for(PermDevice);
}

// template instantiations
template void vector_permute<double>(int n, cuda::const_gpu_ptr<double> x, int incx,
				     cuda::gpu_ptr<double> y, int incy,
				     int const* Perm);

// vector_parallel

template <typename T>
__global__ void cuda_vector_parallel(unsigned n, T const* x, unsigned incx,
				     T const* y, unsigned incy,
				     T* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = x[i*incx] * y[i*incy];
}

template <>
__global__ void cuda_vector_parallel<cuFloatComplex>(unsigned n, cuFloatComplex const* x, unsigned incx,
				     cuFloatComplex const* y, unsigned incy,
				     cuFloatComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCmulf(x[i*incx], y[i*incy]);
}

template <>
__global__ void cuda_vector_parallel<cuDoubleComplex>(unsigned n, cuDoubleComplex const* x, unsigned incx,
				     cuDoubleComplex const* y, unsigned incy,
				     cuDoubleComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCmul(x[i*incx], y[i*incy]);
}

template <typename T>
__global__ void cuda_vector_add_parallel(unsigned n, T const* x, unsigned incx,
					 T const* y, unsigned incy,
					 T* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] += x[i*incx] * y[i*incy];
}



template <>
__global__ void cuda_vector_add_parallel<cuFloatComplex>(unsigned n, cuFloatComplex const* x, unsigned incx,
					 cuFloatComplex const* y, unsigned incy,
					 cuFloatComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCfmaf(x[i*incx], y[i*incy], z[i*incz]);
}

template <>
__global__ void cuda_vector_add_parallel<cuDoubleComplex>(unsigned n, cuDoubleComplex const* x, unsigned incx,
					 cuDoubleComplex const* y, unsigned incy,
					 cuDoubleComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCfma(x[i*incx], y[i*incy], z[i*incz]);
}

template <typename T>
__global__ void cuda_vector_parallel_scaled(unsigned n, T alpha, T const* x, unsigned incx,
					    T const* y, unsigned incy,
					    T* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = alpha * x[i*incx] * y[i*incy];
}

template <>
__global__ void cuda_vector_parallel_scaled<cuFloatComplex>(unsigned n, cuFloatComplex alpha,
                                                            cuFloatComplex const* x, unsigned incx,
                                                            cuFloatComplex const* y, unsigned incy,
                                                            cuFloatComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCmulf(cuCmulf(alpha, x[i*incx]), y[i*incy]);
}

template <>
__global__ void cuda_vector_parallel_scaled<cuDoubleComplex>(unsigned n, cuDoubleComplex alpha,
                                                            cuDoubleComplex const* x, unsigned incx,
                                                            cuDoubleComplex const* y, unsigned incy,
                                                            cuDoubleComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCmul(cuCmul(alpha, x[i*incx]), y[i*incy]);
}

template <typename T>
__global__ void cuda_vector_add_parallel_scaled(unsigned n, T alpha, T const* x, unsigned incx,
						T const* y, unsigned incy,
						T* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] += alpha * x[i*incx] * y[i*incy];
}

template <>
__global__ void cuda_vector_add_parallel_scaled<cuFloatComplex>(unsigned n, cuFloatComplex alpha,
                                                                cuFloatComplex const* x, unsigned incx,
                                                                cuFloatComplex const* y, unsigned incy,
                                                                cuFloatComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCfmaf(cuCmulf(alpha, x[i*incx]), y[i*incy], z[i*incz]);
}

template <>
__global__ void cuda_vector_add_parallel_scaled<cuDoubleComplex>(unsigned n, cuDoubleComplex alpha,
                                                                 cuDoubleComplex const* x, unsigned incx,
                                                                 cuDoubleComplex const* y, unsigned incy,
                                                                 cuDoubleComplex* z, unsigned incz)
{
   unsigned i = blockIdx.x*blockDim.x + threadIdx.x;
   if (i < n)
      z[i*incz] = cuCfma(cuCmul(alpha, x[i*incx]), y[i*incy], z[i*incz]);
}

template <typename T>
void
vector_parallel(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::const_gpu_ptr<T> y, int incy,
		cuda::gpu_ptr<T> z, int incz)
{
   unsigned threads_per_block = 64;
   unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
   cuda_vector_parallel<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>
      (n, cuda_complex_cast(x.device_ptr()), incx,
       cuda_complex_cast(y.device_ptr()), incy,
       cuda_complex_cast(z.device_ptr()), incz);
}

template <typename T>
void
vector_add_parallel(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::const_gpu_ptr<T> y, int incy,
		    cuda::gpu_ptr<T> z, int incz)
{
   unsigned threads_per_block = 64;
   unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
   cuda_vector_add_parallel<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>
      (n, cuda_complex_cast(x.device_ptr()), incx,
       cuda_complex_cast(y.device_ptr()), incy,
       cuda_complex_cast(z.device_ptr()), incz);
}

template <typename T>
void
vector_parallel_scaled(int n, T const& alpha, cuda::const_gpu_ptr<T> x, int incx,
		       cuda::const_gpu_ptr<T> y, int incy,
		       cuda::gpu_ptr<T> z, int incz)
{
   unsigned threads_per_block = 64;
   unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
   cuda_vector_parallel_scaled<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>
      (n, *cuda_complex_cast(&alpha), cuda_complex_cast(x.device_ptr()), incx,
       cuda_complex_cast(y.device_ptr()), incy,
       cuda_complex_cast(z.device_ptr()), incz);
}

template <typename T>
void
vector_add_parallel_scaled(int n, T const& alpha, cuda::const_gpu_ptr<T> x, int incx,
			   cuda::const_gpu_ptr<T> y, int incy,
			   cuda::gpu_ptr<T> z, int incz)
{
   unsigned threads_per_block = 64;
   unsigned nblocks = (n + threads_per_block-1) / threads_per_block;
   cuda_vector_add_parallel_scaled<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>
      (n, *cuda_complex_cast(&alpha), cuda_complex_cast(x.device_ptr()), incx,
       cuda_complex_cast(y.device_ptr()), incy,
       cuda_complex_cast(z.device_ptr()), incz);
}


// template instantiations

template void vector_parallel<float>(int n, cuda::const_gpu_ptr<float> x, int incx,
                                      cuda::const_gpu_ptr<float> y, int incy,
                                      cuda::gpu_ptr<float> z, int incz);
template void vector_parallel<double>(int n, cuda::const_gpu_ptr<double> x, int incx,
                                      cuda::const_gpu_ptr<double> y, int incy,
                                      cuda::gpu_ptr<double> z, int incz);
template void vector_parallel<std::complex<float>>(int n, cuda::const_gpu_ptr<std::complex<float>> x, int incx,
                                                   cuda::const_gpu_ptr<std::complex<float>> y, int incy,
                                                   cuda::gpu_ptr<std::complex<float>> z, int incz);
template void vector_parallel<std::complex<double>>(int n, cuda::const_gpu_ptr<std::complex<double>> x, int incx,
                                                    cuda::const_gpu_ptr<std::complex<double>> y, int incy,
                                                    cuda::gpu_ptr<std::complex<double>> z, int incz);

template void vector_parallel_scaled<float>(int n, float const& alpha, cuda::const_gpu_ptr<float> x, int incx,
                                            cuda::const_gpu_ptr<float> y, int incy,
                                            cuda::gpu_ptr<float> z, int incz);
template void vector_parallel_scaled<double>(int n, double const& alpha, cuda::const_gpu_ptr<double> x, int incx,
                                             cuda::const_gpu_ptr<double> y, int incy,
                                             cuda::gpu_ptr<double> z, int incz);
template void vector_parallel_scaled<std::complex<float>>(int n, std::complex<float> const& alpha,
                                                          cuda::const_gpu_ptr<std::complex<float>> x, int incx,
                                                          cuda::const_gpu_ptr<std::complex<float>> y, int incy,
                                                          cuda::gpu_ptr<std::complex<float>> z, int incz);
template void vector_parallel_scaled<std::complex<double>>(int n, std::complex<double> const& alpha,
                                                           cuda::const_gpu_ptr<std::complex<double>> x, int incx,
                                                           cuda::const_gpu_ptr<std::complex<double>> y, int incy,
                                                           cuda::gpu_ptr<std::complex<double>> z, int incz);

template void vector_add_parallel<float>(int n, cuda::const_gpu_ptr<float> x, int incx,
                                         cuda::const_gpu_ptr<float> y, int incy,
                                         cuda::gpu_ptr<float> z, int incz);
template void vector_add_parallel<double>(int n, cuda::const_gpu_ptr<double> x, int incx,
                                          cuda::const_gpu_ptr<double> y, int incy,
                                          cuda::gpu_ptr<double> z, int incz);
template void vector_add_parallel<std::complex<float>>(int n, cuda::const_gpu_ptr<std::complex<float>> x, int incx,
                                                       cuda::const_gpu_ptr<std::complex<float>> y, int incy,
                                                       cuda::gpu_ptr<std::complex<float>> z, int incz);
template void vector_add_parallel<std::complex<double>>(int n, cuda::const_gpu_ptr<std::complex<double>> x, int incx,
                                                        cuda::const_gpu_ptr<std::complex<double>> y, int incy,
                                                        cuda::gpu_ptr<std::complex<double>> z, int incz);

template void vector_add_parallel_scaled<float>(int n, float const& alpha, cuda::const_gpu_ptr<float> x, int incx,
                                                cuda::const_gpu_ptr<float> y, int incy,
                                                cuda::gpu_ptr<float> z, int incz);
template void vector_add_parallel_scaled<double>(int n, double const& alpha, cuda::const_gpu_ptr<double> x, int incx,
                                          cuda::const_gpu_ptr<double> y, int incy,
                                          cuda::gpu_ptr<double> z, int incz);
template void vector_add_parallel_scaled<std::complex<float>>(int n, std::complex<float> const& alpha,
                                                              cuda::const_gpu_ptr<std::complex<float>> x, int incx,
                                                              cuda::const_gpu_ptr<std::complex<float>> y, int incy,
                                                              cuda::gpu_ptr<std::complex<float>> z, int incz);
template void vector_add_parallel_scaled<std::complex<double>>(int n, std::complex<double> const& alpha,
                                                               cuda::const_gpu_ptr<std::complex<double>> x, int incx,
                                                               cuda::const_gpu_ptr<std::complex<double>> y, int incy,
                                                               cuda::gpu_ptr<std::complex<double>> z, int incz);


__global__ void cuda_dgmm(int M, int N, double const* x, int incx,
                          cuDoubleComplex const* B, int ldb,
                          cuDoubleComplex* C, int ldc)
{
   unsigned Row = blockIdx.x*blockDim.x + threadIdx.x;
   unsigned Col = blockIdx.y*blockDim.y + threadIdx.y;
   if (Row < M && Col < N)
   {
      C[ldc*Col + Row].x = x[Row*incx] * B[ldc*Col + Row].x;
      C[ldc*Col + Row].y = x[Row*incx] * B[ldc*Col + Row].y;
   }
}

void
dgmm(int M, int N,
     cuda::const_gpu_ptr<double> x, int incx,
     cuda::const_gpu_ptr<std::complex<double>> B, int ldb,
     cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   C.wait_for(x);
   C.wait_for(B);
   TRACE_CUDA("cuda_dgmm")(M)(N)(x.device_ptr())(incx)(B.device_ptr())(ldb)(C.device_ptr())(ldc);
   dim3 BlockDim(16, 16);
   dim3 GridDim((M+BlockDim.x-1) / BlockDim.x, (N+BlockDim.y-1) / BlockDim.y);
   cuda_dgmm<<<GridDim, BlockDim, 0, C.get_stream().raw_stream()>>>(M, N, x.device_ptr(), incx,
                                                                    cuda_complex_cast(B.device_ptr()), ldb,
                                                                    cuda_complex_cast(C.device_ptr()), ldc);
   x.wait_for(C);
   B.wait_for(C);
}

} // namespace cuda
