// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/vecptriterator.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
/* -*- C++ -*- $Id$

  vecptriterator.h

  a linear-algebra iterator for a pointer.  This can have either fixed or variable
  stride, specified by the second template parameter which defaults to fixed stride 1.

  To specialize functions for stride 1, use VecPtrIterator<Scalar>.
  To specialize functions for compile-time but fixed stride, use
  VecPtrIterator<Scalar, boost::mpl::int_<int N> >.
  To specialize for variable stride, use VecPtrIterator<Scalar, tagVariable>.

  TODO: is it worthwhile adding fixed size variants too?
*/

#if !defined(VECTORITERATOR_H_JKDSHCORIYT98Y43)
#define VECTORITERATOR_H_JKDSHCORIYT98Y43

#include "interface.h"
#include <boost/type_traits.hpp>
#include <boost/mpl/int.hpp>

namespace LinearAlgebra
{

// A size or stride can either be fixed, in which case
// the type is boost::mpl::int_, or it can be variable
// in which case the type is a dummy tag.
struct tagVariable {};

template <typename Scalar, typename Stride = boost::mpl::int_<1> >
class VecPtrIterator
{
   private:
      static int const Stride_ = Stride::value;
   public:
      // FIXME: should this be make_value ?  Does it matter?
      typedef typename boost::remove_const<Scalar>::type value_type;
      typedef Scalar& reference;
      typedef Scalar* pointer;
      typedef vector_iterator_dense category;

      VecPtrIterator() {}

   VecPtrIterator(Scalar* Base, size_type Size, size_type Index) //__attribute__((always_inline))
         : Base_(Base), Size_(Size), Index_(Index)
     { }

      VecPtrIterator(Scalar* Base, size_type Size, size_type Index,
                     difference_type Str)
         : Base_(Base), Size_(Size), Index_(Index)
     { DEBUG_PRECONDITION_EQUAL(Str, Stride::value); }

      // this overload handles const conversions
      VecPtrIterator(VecPtrIterator<value_type, Stride> const& i)
         : Base_(i.Base_), Size_(i.Size_), Index_(i.Index_) {}

      VecPtrIterator& operator++()
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); ++Index_; return *this; }

      VecPtrIterator operator++(int)
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_);
           return VecPtrIterator(Base_, Size_, Index_++); }

      VecPtrIterator& operator--()
         { DEBUG_CHECK(Index_ > 0)(Index_); --Index_; return *this; }

      VecPtrIterator operator--(int)
         { DEBUG_CHECK(Index_ > 0)(Index_);
           return VecPtrIterator(Base_, Size_, Index_--); }

      VecPtrIterator& operator+=(difference_type n)
         { DEBUG_CHECK(difference_type(Index_) + n >= 0 &&
                       difference_type(Index_) + n <= Size_)(Index_)(Size_);
         Index_ += n; return *this; }

      VecPtrIterator& operator-=(difference_type n)
         { DEBUG_CHECK(difference_type(Index_) - n >= 0 &&
                       difference_type(Index_) - n <= Size_)(Index_)(Size_);
         Index_ -= n; return *this; }

      size_type index() const { return Index_; }

      reference operator*() const
      { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); return Base_[Index_ * Stride_]; }

      pointer operator->() const
      {  DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); return Base_ + Index_ * Stride_; }

      reference operator[](difference_type n) const
      {  DEBUG_RANGE_CHECK_OPEN(difference_type(Index_) + n, 0, difference_type(Size_));
         return Base_[(difference_type(Index_)+n) * Stride_]; }

      operator bool() const { return Index_ != Size_; }

      Scalar* ptr() const { return Base_ + difference_type(Index_) * Stride_; }
      Scalar* end() const { return Base_ + difference_type(Size_) * Stride_; }
      size_type size() const { return Size_ - Index_; }
      difference_type stride() const { return Stride_; }

      static difference_type const static_stride = Stride_;

      // For BLAS, the required pointer is to the lowest element in memory.
      // For negative stride,
      // this is actually the 'back' of the array.
      Scalar* blas_base() const
      { return Base_ + difference_type(Stride_ > 0 ? Index_ : Size_) * Stride_; }

   private:
      Scalar* Base_;
      size_type Size_, Index_;

   template <typename Other, typename S2> friend class VecPtrIterator;
};

template <typename Scalar>
class VecPtrIterator<Scalar, tagVariable>
{
   public:
      typedef typename boost::remove_const<Scalar>::type value_type;
      typedef Scalar& reference;
      typedef Scalar* pointer;
      typedef vector_iterator_dense category;

      VecPtrIterator() {}

      VecPtrIterator(VecPtrIterator<value_type, tagVariable> const& i)
         : Base_(i.Base_), Size_(i.Size_), Index_(i.Index_), Stride_(i.Stride_) {}

      VecPtrIterator(Scalar* Base, size_type Size, size_type Index, difference_type Stride)
         : Base_(Base), Size_(Size), Index_(Index), Stride_(Stride) { }

      VecPtrIterator& operator++()
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); ++Index_; return *this; }
      //         { ++Index_; return *this; }

      VecPtrIterator operator++(int)
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_);
           return VecPtrIterator(Base_, Size_, Index_++); }

      VecPtrIterator& operator--()
         { DEBUG_CHECK(Index_ > 0)(Index_); --Index_; return *this; }

      VecPtrIterator operator--(int)
         { DEBUG_CHECK(Index_ > 0)(Index_);
           Index_ -= Stride_; return VecPtrIterator(Base_, Size_, Index_+Stride_); }

      VecPtrIterator& operator+=(difference_type n)
         { DEBUG_CHECK(difference_type(Index_) + n >= 0 &&
                       difference_type(Index_) + n <= Size_)(Index_)(Size_);
         Index_ += n; return *this; }

      VecPtrIterator& operator-=(difference_type n)
         { DEBUG_CHECK(difference_type(Index_) - n >= 0 &&
                       difference_type(Index_) - n <= Size_)(Index_)(Size_);
         Index_ -= n; return *this; }

      size_type index() const { return Index_; }

      reference operator*() const
      { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); return Base_[Index_ * Stride_]; }
   //      { return Base_[Index_ * Stride_]; }

      pointer operator->() const
      {  DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); return Base_ + Index_ * Stride_; }

      reference operator[](difference_type n) const
      {  DEBUG_RANGE_CHECK_OPEN(difference_type(Index_) + n, 0, difference_type(Size_));
         return Base_[(difference_type(Index_)+n) * Stride_]; }

      operator bool() const { return Index_ != Size_; }

      Scalar* ptr() const { return Base_ + difference_type(Index_) * Stride_; }
      Scalar* end() const { return Base_ + difference_type(Size_) * Stride_; }
      size_type size() const { return Size_ - Index_; }
      difference_type stride() const { return Stride_; }

      // For BLAS, the required pointer is to the lowest element in memory.
      // For negative stride,
      // this is actually the 'back' of the array.
      Scalar* blas_base() const
      { return Base_ + difference_type(Stride_ > 0 ? Index_ : Size_) * Stride_; }

   private:
      Scalar* Base_;
      size_type Size_, Index_;
      difference_type Stride_;

   template <typename Other, typename S2> friend class VecPtrIterator;
};

// this specialization is probably overly-paranoid
template <typename Scalar>
class VecPtrIterator<Scalar, boost::mpl::int_<1> >
{
   private:
      static int const Stride_ = 1;
   public:
      typedef typename boost::remove_const<Scalar>::type value_type;
      typedef Scalar& reference;
      typedef Scalar* pointer;
      typedef vector_iterator_dense category;

      VecPtrIterator() {}

      // this overload handles const conversions
      VecPtrIterator(VecPtrIterator<value_type, boost::mpl::int_<1> > const& i)
         : Base_(i.Base_), Size_(i.Size_), Index_(i.Index_) {}

      VecPtrIterator(Scalar* Base, size_type Size, size_type Index,
                     difference_type Stride = 1)
         : Base_(Base), Size_(Size), Index_(Index) { DEBUG_PRECONDITION_EQUAL(Stride, 1); }

      VecPtrIterator& operator++()
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); ++Index_; return *this; }

      VecPtrIterator operator++(int)
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_);
           return VecPtrIterator(Base_, Size_, Index_++); }

      VecPtrIterator& operator--()
         { DEBUG_CHECK(Index_ > 0)(Index_); --Index_; return *this; }

      VecPtrIterator operator--(int)
         { DEBUG_CHECK(Index_ > 0)(Index_);
           Index_ -= Stride_; return VecPtrIterator(Base_, Size_, Index_+Stride_); }

      VecPtrIterator& operator+=(difference_type n)
         { DEBUG_CHECK(difference_type(Index_) + n >= 0 &&
                       Index_ + size_type(n) <= Size_)(Index_)(Size_);
         Index_ += n; return *this; }

      VecPtrIterator& operator-=(difference_type n)
         { DEBUG_CHECK(difference_type(Index_) - n >= 0 &&
                       difference_type(Index_) - n <= Size_)(Index_)(Size_);
         Index_ -= n; return *this; }

      size_type index() const { return Index_; }

      reference operator*() const
      { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); return Base_[Index_]; }

      pointer operator->() const
      {  DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); return Base_ + Index_; }

      reference operator[](difference_type n) const
      {  DEBUG_RANGE_CHECK_OPEN(difference_type(Index_) + n, 0, difference_type(Size_));
         return Base_[difference_type(Index_)+n]; }

      operator bool() const { return Index_ != Size_; }

      Scalar* ptr() const { return Base_ + Index_; }
      Scalar* end() const { return Base_ + Size_; }
      size_type size() const { return Size_ - Index_; }
      difference_type stride() const { return 1; }

      static difference_type const static_stride = 1;

      Scalar* blas_base() const { return Base_ + Index_; }

   private:
      Scalar* Base_;
      size_type Size_, Index_;

   template <typename Other, typename S2> friend class VecPtrIterator;
};

} // namespace LinearAlgebra

#endif
