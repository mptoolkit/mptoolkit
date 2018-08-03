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

#if !defined(MPTOOLKIT_BLAS_RANGE_H)
#define MPTOOLKIT_BLAS_RANGE_H

#include "vectorref.h"

namespace blas
{

//
// Range
//
// A Range represents some segment of another vector.
// It is itself a dense vector, as well as a unary function.
//

class Range;
class cpu_tag;

class RangeIterator
{
   public:
      using value_type      = int;
      using reference       = value_type const&;
      using pointer         = value_type const*;
      using category        = std::random_access_iterator_tag;

      RangeIterator() {}

      RangeIterator(int n_)
         : n(n_) {}

      RangeIterator& operator++() { ++n; return *this; }
      RangeIterator operator++(int) { return RangeIterator(n++); }
      RangeIterator& operator+=(int i) { n += i; return *this; }

      value_type operator[](int i) const { return n+i; }

      value_type operator*() const { return n; }

      pointer operator*() { return &n; }

   private:
      int n;
};

class Range : public VectorRef<int, Range, cpu_tag>
{
   public:
      using value_type     = int;
      using const_iterator = RangeIterator;
      using iterator       = RangeIterator;

      Range() : First_(0), Last_(0) {}

      Range(int First, int Last)
         : First_(First), Last_(Last) { DEBUG_PRECONDITION(First <= Last); }

      Range(Range&& Other) = default;

      Range(Range const& Other) : VectorRef<int, Range, cpu_tag>(),
                                  First_(Other.First_), Last_(Other.Last_) {}

      Range& operator=(Range const& r) { First_ = r.First_; Last_ = r.Last_; return *this; }

      const_iterator begin() const { return iterator(First_); }
      const_iterator cbegin() const { return iterator(First_); }

      const_iterator end() const { return iterator(First_); }
      const_iterator cend() const { return iterator(First_); }

      int first() const { return First_; }
      int last() const { return Last_; }

      int size() const { return Last_-First_; }

      int operator()(int n) const { return n + First_; }
      int operator[](int n) const { return n + First_; }

      static constexpr int stride() { return 1; }

      bool operator==(Range const& r) const { return First_ == r.First_ && Last_ == r.Last_; }
      bool operator!=(Range const& r) const { return First_ != r.First_ || Last_ != r.Last_; }

   private:
      int First_, Last_;
};

inline
Range range(int first, int last)
{
   return Range(first, last);
}

struct all_t {};

extern all_t all;

inline
std::ostream& operator<<(std::ostream& out, Range const& r)
{
   out << "range[" << r.first() << "," << r.last() << ']';
   return out;
}

} // namespace blas

#endif
