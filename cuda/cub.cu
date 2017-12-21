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
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx, cuda::gpu_ref<T>&& r)
{
   r.wait(x.sync());
   std::size_t TempStorageBytes = 0;
   void* TempStorage = nullptr;
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.get_stream().raw_stream());
   TempStorage = cuda::allocate_gpu_temporary(TempStorageBytes);
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.get_stream().raw_stream());
   cuda::free_gpu_temporary(TempStorage, TempStorageBytes);
   x.wait(r.sync());
}

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx, cuda::gpu_ref<T>& r)
{
   r.wait(x.sync());
   std::size_t TempStorageBytes = 0;
   void* TempStorage = nullptr;
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.get_stream().raw_stream());
   TempStorage = cuda::allocate_gpu_temporary(TempStorageBytes);
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, detail::stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.get_stream().raw_stream());
   cuda::free_gpu_temporary(TempStorage, TempStorageBytes);
   x.wait(r.sync());
}

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx, cuda::gpu_ref<std::complex<T>>&& r)
{
   cuda::stream Stream1, Stream2;
   std::size_t TempStorageBytes1 = 0;
   void* TempStorage1 = nullptr;
   // real part
   // workspace query
   Stream1.wait(x.sync());
   Stream1.wait(r.sync());
   cub::DeviceReduce::Sum(TempStorage1, TempStorageBytes1, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr())), Size, r.get_stream().raw_stream());
   TempStorage1 = cuda::allocate_gpu_temporary(TempStorageBytes1);
   // do the work
   cub::DeviceReduce::Sum(TempStorage1, TempStorageBytes1, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr())), Size, r.get_stream().raw_stream());

   std::size_t TempStorageBytes2 = 0;
   void* TempStorage2 = nullptr;

   // imag part
   // workspace query
   Stream2.wait(x.sync());
   Stream2.wait(r.sync());
   cub::DeviceReduce::Sum(TempStorage2, TempStorageBytes2, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size, r.get_stream().raw_stream());
   TempStorage1 = cuda::allocate_gpu_temporary(TempStorageBytes1);
   // do the work
   cub::DeviceReduce::Sum(TempStorage2, TempStorageBytes2, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size, r.get_stream().raw_stream());

   r.wait(Stream1.record());
   r.wait(Stream2.record());

   cuda::free_gpu_temporary(TempStorage1, TempStorageBytes1);
   cuda::free_gpu_temporary(TempStorage2, TempStorageBytes2);
}

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<std::complex<T>> const& x, int incx, cuda::gpu_ref<std::complex<T>>& r)
{
   cuda::stream Stream1, Stream2;
   std::size_t TempStorageBytes1 = 0;
   void* TempStorage1 = nullptr;
   // real part
   // workspace query
   Stream1.wait(x.sync());
   Stream1.wait(r.sync());
   cub::DeviceReduce::Sum(TempStorage1, TempStorageBytes1, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr())), Size, r.get_stream().raw_stream());
   TempStorage1 = cuda::allocate_gpu_temporary(TempStorageBytes1);
   // do the work
   cub::DeviceReduce::Sum(TempStorage1, TempStorageBytes1, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr())), incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr())), Size, r.get_stream().raw_stream());

   std::size_t TempStorageBytes2 = 0;
   void* TempStorage2 = nullptr;

   // imag part
   // workspace query
   Stream2.wait(x.sync());
   Stream2.wait(r.sync());
   cub::DeviceReduce::Sum(TempStorage2, TempStorageBytes2, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size, r.get_stream().raw_stream());
   TempStorage1 = cuda::allocate_gpu_temporary(TempStorageBytes1);
   // do the work
   cub::DeviceReduce::Sum(TempStorage2, TempStorageBytes2, 
			  detail::stride_ptr<T const>(static_cast<T const*>(static_cast<void const*>(x.device_ptr()))+1, incx*2),
                          static_cast<T*>(static_cast<void*>(r.device_ptr()))+1, Size, r.get_stream().raw_stream());

   r.wait(Stream1.record());
   r.wait(Stream2.record());

   cuda::free_gpu_temporary(TempStorage1, TempStorageBytes1);
   cuda::free_gpu_temporary(TempStorage2, TempStorageBytes2);
}

// template instantiations
template void vector_sum<double>(int Size, cuda::const_gpu_ptr<double> const& x, int incx, cuda::gpu_ref<double>&& r);
template void vector_sum<double>(int Size, cuda::const_gpu_ptr<double> const& x, int incx, cuda::gpu_ref<double>& r);

template void vector_sum<double>(int Size, 
				 cuda::const_gpu_ptr<std::complex<double>> const& x, int incx, 
				 cuda::gpu_ref<std::complex<double>>&& r);
template void vector_sum<double>(int Size, 
				 cuda::const_gpu_ptr<std::complex<double>> const& x, int incx, 
				 cuda::gpu_ref<std::complex<double>>& r);

template void vector_sum<float>(int Size, cuda::const_gpu_ptr<float> const& x, int incx, cuda::gpu_ref<float>&& r);
template void vector_sum<float>(int Size, cuda::const_gpu_ptr<float> const& x, int incx, cuda::gpu_ref<float>& r);

template void vector_sum<float>(int Size, 
				cuda::const_gpu_ptr<std::complex<float>> const& x, int incx, 
				cuda::gpu_ref<std::complex<float>>&& r);
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
   int* PermDevice = static_cast<int*>(allocate_gpu_temporary(n*sizeof(int)));
   memcpy_host_to_device_async(y.get_stream(), Perm, PermDevice, n);
   if (incx == 1 && incy == 1)
   {
      unsigned threads_per_block = 64;
      unsigned nblocks = (n + threads_per_block-1) / threads_per_block;      
      permute_vector_cuda<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>(n, x.device_ptr(), y.device_ptr(), 
											 PermDevice);
   }
   else
   {
      unsigned threads_per_block = 64;
      unsigned nblocks = (n + threads_per_block-1) / threads_per_block;      
      permute_vector_stride_cuda<<<nblocks,threads_per_block, 0, y.get_stream().raw_stream()>>>(n, x.device_ptr(), incx,
												y.device_ptr(), incy,
												PermDevice);
   }
   x.wait_for(y);
   free_gpu_temporary(static_cast<void*>(PermDevice), n*sizeof(int));
}

// template instantiations
template void vector_permute<double>(int n, cuda::const_gpu_ptr<double> x, int incx, 
				     cuda::gpu_ptr<double> y, int incy, 
				     int const* Perm);

} // namespace cuda
