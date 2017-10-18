// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/slice.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_SLICE_H)
#define MPTOOLKIT_BLAS_SLICE_H

#include "vector.h"

namespace blas
{


//
// Slice
//
// A Slice represents a section view of a container, slice(i) = start + stride * i.
// It is itself a vector, implementing the VectorConstStride interface.
//

class SliceIterator
{
   public:
      using value_type      = int;
      using reference       = value_type const&;
      using pointer         = value_type const*;
      using category        = std::random_access_iterator_tag;

      SliceIterator() {}
      SliceIterator(int Stride, int n = 0)
         : Start_(Start), Size_(Size), Stride_(Stride), n_(n) {}

      SliceIterator& operator++() { n_ += Stride_; return *this; }
      SliceIterator operator++(int) { int n_old = n_; n_ += Stride; return SliceIterator(Stride_, n_); }

      SliceIterator& operator--() { n_ -= Stride_; return *this; }
      SliceIterator operator--(int) { int n_old = n_; n_ -= Stride; return SliceIterator(Stride_, n_); }

      SliceIterator& operator+=(int i) { n_ += i*Stride_; return *this; }
      SliceIterator& operator-=(int i) { n_ -= i*Stride_; return *this; }

      int operator*() const { return n_; }

      int const* operator->() const { return &n_; }

      int operator[](int n) const
      {
         return n_ + n*Stride_;
      }

      int stride() const { return Stride_; }

   private:
      int Stride_;
      int n_;
};

class Slice : public VectorRef<int, Slice, cpu_tag>
{
   public:
      using value_type     = int;
      using const_iterator = SliceIterator;
      using iterator       = SliceIterator;

      Slice() : Start_(0), Size_(0), Stride_(0) {}
      Slice(int Start, int Size, int Stride)
        : Start_(Start), Size_(Size), Stride_(Stride) {}

      int start() const { return Start_; }
      int size() const { return Size_; }
      int stride() const { return Stride_; }

      int operator()(int n) const
      { return int(Start_) + int(n)*Stride_; }

      int operator[](int n) const
      { return int(Start_) + int(n)*Stride_; }

      bool operator==(Slice const& s) const
      { return Start_ == s.Start_ && Size_ == s.Size_ && Stride_ == s.Stride_; }

      bool operator!=(Slice const& s) const
      { return Start_ != s.Start_ || Size_ != s.Size_ || Stride_ != s.Stride_; }

      const_iterator begin() const
      { return const_iterator(Start_, Stride_); }
      const_iterator cbegin() const
      { return const_iterator(Start_, Stride_); }

      const_iterator end() const
      { return const_iterator(Start_ + Size_*Stride_, Stride_); }
      const_iterator cend() const
      { return const_iterator(Start_ + Size_*Stride_, Stride_); }

   private:
      int Start_, Size_;
      int Stride_;
};

inline
Slice slice(int start, int size, int stride)
{
   return Slice(start, size, stride);
}


} // namespace blas

#endif
