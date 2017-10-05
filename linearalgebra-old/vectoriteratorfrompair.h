// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectoriteratorfrompair.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  vectoriteratorfrompair.h

  An adaptor to construct a linear algebra sparse iterator from
  a standard iterator over pair<size_type, T>

  Created 2005-01-06 Ian McCulloch
*/

#if !defined(VECTORITERATORFROMPAIR_H_JKCHWUITY438T97UY89PU908UWP)
#define VECTORITERATORFROMPAIR_H_JKCHWUITY438T97UY89PU908UWP

#include <iterator>
#include "common/trace.h"
#include "interface.h"

namespace LinearAlgebra
{

// VectorIteratorFromPair is a generic sparse iterator from a standard
// iterator over std::pair<size_type, T> or similar.
template <typename BaseIter, typename Category>
class VectorConstIteratorFromPair;

template <typename BaseIter, typename Category>
class VectorIteratorFromPair
{
   private:
      typedef BaseIter base_iter_type;
      typedef typename std::iterator_traits<base_iter_type>::value_type base_value_type;

   public:
      typedef typename base_value_type::second_type value_type;
      typedef value_type& reference;
      typedef value_type* pointer;
      typedef Category category;

      VectorIteratorFromPair() {}

      VectorIteratorFromPair(BaseIter const& Here, BaseIter const& End)
         : Here_(Here), End_(End) {}

      VectorIteratorFromPair& operator++()
      { DEBUG_CHECK(Here_ != End_); ++Here_; return *this; }

      VectorIteratorFromPair operator++(int)
      { DEBUG_CHECK(Here_ != End_); return VectorIteratorFromPair(Here_++, End_); }

      // we cannot do a debug check on operator-- without storing the begin
      // iterator as well; but we don't need begin for any other purpose?
      VectorIteratorFromPair& operator--()
      { --Here_; return *this; }

      VectorIteratorFromPair operator--(int)
      { return VectorIteratorFromPair(Here_--, End_); }

      size_type index() const { DEBUG_CHECK(Here_ != End_); return Here_->first; }

      reference operator*() const { return Here_->second; }

      pointer operator->() const { return &Here_->second; }

      operator bool() const { return Here_ != End_; }

   private:
      base_iter_type Here_, End_;

   template <typename OtherBaseIter, typename OtherCat>
   friend class VectorConstIteratorFromPair;
};

template <typename BaseIter, typename Category>
class VectorConstIteratorFromPair
{
   private:
      typedef BaseIter base_iter_type;
      typedef typename std::iterator_traits<base_iter_type>::value_type base_value_type;

   public:
      typedef typename base_value_type::second_type value_type;
      typedef value_type const& reference;
      typedef value_type const* pointer;
      typedef Category category;

      VectorConstIteratorFromPair() {}

      VectorConstIteratorFromPair(BaseIter const& Here, BaseIter const& End)
         : Here_(Here), End_(End) {}

      template <typename OtherBaseIter>
      VectorConstIteratorFromPair(VectorIteratorFromPair
                                      <OtherBaseIter, Category> const& Other)
         : Here_(Other.Here_), End_(Other.End_) {}

      VectorConstIteratorFromPair& operator++()
      { DEBUG_CHECK(Here_ != End_); ++Here_; return *this; }

      VectorConstIteratorFromPair operator++(int)
      { DEBUG_CHECK(Here_ != End_); return VectorConstIteratorFromPair(Here_++, End_); }

      // we cannot do a debug check on operator-- without storing the begin
      // iterator as well; but we don't need begin for any other purpose?
      VectorConstIteratorFromPair& operator--()
      { --Here_; return *this; }

      VectorConstIteratorFromPair operator--(int)
      { return VectorConstIteratorFromPair(Here_--, End_); }

      size_type index() const { DEBUG_CHECK(Here_ != End_); return Here_->first; }

      reference operator*() const { return Here_->second; }

      pointer operator->() const { return &Here_->second; }

      operator bool() const { return Here_ != End_; }

   private:
      base_iter_type Here_, End_;
};

} // namespace LinearAlgebra

#endif
