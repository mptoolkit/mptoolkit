// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/detail/sparserowiterator.h
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

#if !defined(MPTOOLKIT_BLAS_DETAIL_SPARSEROWERATOR_H)
#define MPTOOLKIT_BLAS_DETAIL_SPARSEROWERATOR_H

#include "common/trace.h"
#include <map>

namespace blas
{

namespace detail
{

// proxy class for an element of a sparse row vector (indexed by a column)
template <typename T>
struct SparseElementRowRef
{
   int index;
   T&  value;

   SparseElementRowRef() = delete;

   SparseElementRowRef(int Col_, T& Value_) : index(Col_), value(Value_) {}

   int col() const { return index; }

   operator T&() const { return value; }

};

template <typename T>
class SparseRowIterator
{
   public:
      using base_iterator = typename std::map<int, T>::iterator;

      using value_type = SparseElementRowRef<T>;
      using pointer    = T*;
      using reference  = SparseElementRowRef<T>;

      using iterator_category = typename base_iterator::iterator_category;
      using difference_type   = typename base_iterator::difference_type;

      SparseRowIterator() = default;

      explicit SparseRowIterator(base_iterator Other) : I(std::move(Other)) {}

      SparseRowIterator(SparseRowIterator const& Other) = default;
      SparseRowIterator(SparseRowIterator&& Other) = default;

      SparseRowIterator& operator=(SparseRowIterator const& Other) = default;
      SparseRowIterator& operator=(SparseRowIterator&& Other) = default;

      int col() const { return I->first; }

      int index() const { return I->first; }

      T& value() const { return I->second; }

      reference operator*() const { return reference(I->first, I->second); }
      pointer operator->() const { return &I->second; }

      SparseRowIterator& operator++() noexcept
      {
         ++I;
         return *this;
      }

      SparseRowIterator operator++(int) noexcept
      {
         return SparseRowIterator(I++);
      }

      SparseRowIterator& operator--() noexcept
      {
         --I;
         return *this;
      }

      SparseRowIterator operator--(int) noexcept
      {
         return SparseRowIterator(I--);
      }

      bool operator==(SparseRowIterator const& x) const
      {
         return I == x.I;
      }

      bool operator!=(SparseRowIterator const& x) const
      {
         return I != x.I;
      }

      base_iterator const& base() const { return I; }
      base_iterator& base() { return I; }

   private:
      base_iterator I;
};

template <typename T>
class ConstSparseRowIterator
{
   public:
      using base_iterator = typename std::map<int, T>::const_iterator;

      using value_type = SparseElementRowRef<T const>;
      using pointer    = T const*;
      using reference  = SparseElementRowRef<T const>;

      using iterator_category = typename base_iterator::iterator_category;
      using difference_type   = typename base_iterator::difference_type;

      ConstSparseRowIterator() = default;

      explicit ConstSparseRowIterator(base_iterator Other) : I(std::move(Other)) {}

      ConstSparseRowIterator(ConstSparseRowIterator const& Other) = default;
      ConstSparseRowIterator(ConstSparseRowIterator&& Other) = default;

      ConstSparseRowIterator& operator=(ConstSparseRowIterator const& Other) = default;
      ConstSparseRowIterator& operator=(ConstSparseRowIterator&& Other) = default;

#if 0
      ConstSparseRowIterator(ConstSparseRowIterator<T> const& Other)
         : I(Other.base()) {}

      ConstSparseRowIterator(ConstSparseRowIterator<T>&& Other)
         : I(std::move(Other.base())) {}
#endif

      int col() const { return I->first; }

      int index() const { return I->first; }

      T const& value() const { return I->second; }

      reference operator*() const { return reference(I->first, I->second); }
      pointer operator->() const { return &I->second; }

      ConstSparseRowIterator& operator++() noexcept
      {
         ++I;
         return *this;
      }

      ConstSparseRowIterator operator++(int) noexcept
      {
         return ConstSparseRowIterator(I++);
      }

      ConstSparseRowIterator& operator--() noexcept
      {
         --I;
         return *this;
      }

      ConstSparseRowIterator operator--(int) noexcept
      {
         return ConstSparseRowIterator(I--);
      }

      bool operator==(ConstSparseRowIterator const& x) const
      {
         return I == x.I;
      }

      bool operator!=(ConstSparseRowIterator const& x) const
      {
         return I != x.I;
      }

      base_iterator const& base() const { return I; }
      base_iterator& base() { return I; }

   private:
      base_iterator I;
};

template <typename T>
inline
bool operator==(SparseRowIterator<T> const& x, ConstSparseRowIterator<T> const& y)
{
   return x.base() == y.base();
}

template <typename T>
inline
bool operator==(ConstSparseRowIterator<T> const& x, SparseRowIterator<T> const& y)
{
   return x.base() == y.base();
}

template <typename T>
inline
bool operator!=(SparseRowIterator<T> const& x, ConstSparseRowIterator<T> const& y)
{
   return x.base() != y.base();
}

template <typename T>
inline
bool operator!=(ConstSparseRowIterator<T> const& x, SparseRowIterator<T> const& y)
{
   return x.base() != y.base();
}

} // namespace detail

} // namespace blas

#endif
