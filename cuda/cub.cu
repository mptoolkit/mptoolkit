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

#include "cuda/cub.h"
#include "cub/cub.cuh"
#include <type_traits>

namespace cub
{

template <typename T>
class stride_ptr
{
   public:
      using iterator_category = std::random_access_iterator_tag;
      using value_type        = typename std::remove_const<T>::type;
      using pointer           = T*;
      using reference         = T&;

      stride_ptr() {}
      stride_ptr(T* x_, int stride_) : x(x_), stride(stride_) {}
      stride_ptr(stride_ptr<value_type> const& Other) : x(Other.x), stride(Other.stride) {}

      reference operaor*() const { return *x; }
      pointer operator->() const { return x; }

      reference operator[](int i) const { return x[i*stride]; }

      stride_ptr& operator++() { x += stride; return this; }
      stride_ptr operator++(int) { T* Temp = x; x += stride; return stride_ptr(Temp, stride); }

      stride_ptr& operator--() { x -= stride; return this; }
      stride_ptr operator--(int) { T* Temp = x; x -= stride; return stride_ptr(Temp, stride); }

      stride_ptr& operator+=(int i) { x += i*stride; return *this; }

   private:
      friend class stride_ptr<value_type>;

      T* x;
      int stride;
};

template <typename T>
void
vector_sum(int Size, const_gpu_ptr<T> const& x, int incx, gpu_ref<T>& r)
{
   int TempStorageBytes = 0;
   double* TempStorage = nullptr;
   cub::DeviceReduce::Sum(TempStorage, TempStorageBytes, stride_ptr<T const>(x.device_ptr(), incx),
                          r.device_ptr(), Size, r.stream().raw_stream());
}

// template instantiations
template void vector_sum(int Size, const_gpu_ptr<double> const& x, int incx, gpu_ref<double>& r);

} // namespace cub

