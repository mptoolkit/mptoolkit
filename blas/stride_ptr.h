// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/stride_ptr.h
//
// Copyright (C) 2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// stride_ptr is a pointer type that has a stride;
// it functions like a simple strided iterator over T*
//
// Note that stride iterators in C++ suffer from a dubious status with respect
// to undefined behaviour of the end() iterator.  For a stride > 1, the pointer
// to the end iterator lies at (stride) elements beyond the end of the array.
// The C/C++ standards guarantee that it is possible to form (but not deference!)
// a pointer to one-past-the-end of an array, but not further than that.  In practice,
// POSIX specifies a flat memory model, which gives additional guarantees beyond
// what the C/C++ standards specify.  This is probably still not enough to guarantee
// that a stride_ptr functioning as an end iterator is guaranteed to work, but in
// practice it is probably OK.

#if !defined(MPTOOLKIT_BLAS_STRIDE_PTR_H)
#define MPTOOLKIT_BLAS_STRIDE_PTR_H

#include <iterator>

namespace blas
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

      stride_ptr(T* x, int stride) : x_(x), stride_(stride) {}

      reference operator*() const { return *x_; }
      pointer operator->() const { return x_; }

      reference operator[](int i) const { return x_[i*stride_]; }

      stride_ptr& operator++() { x_ += stride_; return this; }
      stride_ptr operator++(int) { T* Temp = x_; x_ += stride_; return stride_ptr(Temp, stride_); }

      stride_ptr& operator--() { x_ -= stride_; return this; }
      stride_ptr operator--(int) { T* Temp = x_; x_ -= stride_; return stride_ptr(Temp, stride_); }

      stride_ptr& operator+=(int i) { x_ += i*stride_; return *this; }

      int stride() const { return stride_; }

   private:
      T* x_;
      int stride_;
};

} // namespace blas

#endif
