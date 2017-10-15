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

namespace cub
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
      friend class stride_ptr<value_type>;

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

template <typename T>
void
vector_sum(int Size, cuda::const_gpu_ptr<T> const& x, int incx, cuda::gpu_ref<T>& r)
{
   std::size_t TempStorageBytes = 0;
   void* TempStorage = nullptr;
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.get_stream().raw_stream());
   TempStorage = cuda::allocate_gpu_temporary(TempStorageBytes);
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.get_stream().raw_stream());
   cuda::free_gpu_temporary(TempStorage, TempStorageBytes);
}

// template instantiations
template void vector_sum<double>(int Size, cuda::const_gpu_ptr<double> const& x, int incx, cuda::gpu_ref<double>& r);

} // namespace cub

