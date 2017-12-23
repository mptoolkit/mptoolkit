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

      __host__ __device__ stride_ptr() {}
      __host__ __device__ stride_ptr(T* x, int stride) : x_(x), stride_(stride) {}
      __host__ __device__ stride_ptr(stride_ptr const& Other) = default;

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
   std::size_t TempStorageBytes = 0;
   cub::DeviceReduce::Sum(nullptr, TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size);
   gpu_buffer<unsigned char> TempBuffer = allocate_gpu_temporary<unsigned char>(TempStorageBytes);
   TempBuffer.wait_for(x);
   TempBuffer.wait_for(r);
   cub::DeviceReduce::Sum(static_cast<void*>(TempBuffer.device_ptr()), 
			  TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, TempBuffer.get_stream().raw_stream());
   r.wait(TempBuffer.sync());
   x.wait(r.sync());
}

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx, 
	   cuda::gpu_ref<std::complex<T>>& r)
{
   std::size_t TempStorageBytes1 = 0;
   // real part
   // workspace query
   cub::DeviceReduce::Sum(nullptr, TempStorageBytes1, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr())), Size);

   gpu_buffer<unsigned char> TempBuffer1 = allocate_gpu_temporary<unsigned char>(TempStorageBytes1);
   TempBuffer1.wait_for(x);
   TempBuffer1.wait_for(r);
   // do the work
   cub::DeviceReduce::Sum(static_cast<void*>(TempBuffer1.device_ptr()), TempStorageBytes1,
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr())), Size, 
			  TempBuffer1.get_stream().raw_stream());

   // imag part
   // workspace query
   std::size_t TempStorageBytes2 = 0;
   cub::DeviceReduce::Sum(nullptr, TempStorageBytes2, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size);
   gpu_buffer<unsigned char> TempBuffer2 = allocate_gpu_temporary<unsigned char>(TempStorageBytes2);
   TempBuffer2.wait_for(x);
   TempBuffer2.wait_for(r);
   // do the work
   cub::DeviceReduce::Sum(static_cast<void*>(TempBuffer2.device_ptr()), TempStorageBytes2, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size, 
			  TempBuffer2.get_stream().raw_stream());

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

} // namespace cuda
