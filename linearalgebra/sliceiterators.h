// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/sliceiterators.h
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
  sliceiterators.h

  Iterator adaptors for ranges, slices

  Created 2005-01-12 Ian McCulloch

  Still to be done: some tests, and remaining fixed-stride versions.
*/

#if !defined(SLICEITERATORS_H_HCUIRH83479Y98YP)
#define SLICEITERATORS_H_HCUIRH83479Y98YP

namespace LinearAlgebra
{

template <typename BaseIter, typename Stride = boost::mpl::int_<1>,
   typename Category = typename BaseIter::category>
class StrideIterator;

// base case does sparse, injective or hashed.  variable stride.

// FIXME: the rounding for negative strides is probably broken here
template <typename BaseIter, typename Category>
class StrideIterator<BaseIter, tagVariable, Category>
{
   public:
      typedef BaseIter base_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer   pointer;
      typedef typename base_type::category  category;

      StrideIterator(base_type const& Base, size_type Start, size_type Size, 
		     difference_type Stride)
	 : Base_(Base), Start_(Start), Size_(Size), Stride_(Stride) 
	 { while (this->off_stride()) ++Base_; }

      StrideIterator& operator++()
      { DEBUG_CHECK(bool(*this)); ++Base_; while (this->off_stride()) ++Base_; return *this; }

      StrideIterator operator++(int) 
         { DEBUG_CHECK(bool(*this)); 
	   return StrideIterator(Base_++, Start_, Size_, Stride_); }

      // no operator--, can't define it properly

      size_type index() const { return (Base_.index() - Start_) / Stride_; }

      reference operator*() const { return *Base_; }

      pointer operator->() const { return Base_::operator->(); }

      operator bool() const { return Base_; }

      // hash-iterator members	 
      size_type size() const { return std::min(Base_.size(), Size_); }  // hmm
      reference operator()(difference_type n) const
      { return Base_(Start_ + n * Stride_); }
      bool element_exists() const { return Base_.element_exists(Start_ + n * Stride_); }

      base_type const& base() const { return Base_; }
      base_type& base() { return Base_; }
      difference_type stride() const { return Stride_; }

   private:
      bool off_stride() const { return Base_ && (Base_.index() - Start_) % Stride_ != 0; }

      base_type Base_;
      difference_type Start_;
      size_type Size_;
      difference_type Stride_;
};

// we get a speed improvement on the ordered case because we can stop
// once the index gets too big

// FIXME: the rounding for negative strides is probably broken here
template <typename BaseIter>
class StrideIterator<BaseIter, tagVariable, vector_iterator_ordered>
{
   public:
      typedef BaseIter base_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer   pointer;
      typedef typename base_type::category  category;

      StrideIterator(base_type const& Base, size_type Start, size_type Size, difference_type Stride)
	 : Base_(Base), Start_(Start), Size_(Size), Stride_(Stride),
	 LastIndex_(Start_ + Size_ * Stride_)
	 { while (this->off_stride()) ++Base_; }

      StrideIterator& operator++()
      { DEBUG_CHECK(bool(*this)); ++Base_; while (this->off_stride()) ++Base_; return *this; }

      StrideIterator operator++(int) 
         { DEBUG_CHECK(bool(*this)); 
	   return StrideIterator(Base_++, Start_, Size_, Stride_); }

      // no operator--, can't define it properly

      size_type index() const { return (Base_.index() - Start_) / Stride_; }

      reference operator*() const { return *Base_; }

      pointer operator->() const { return Base_::operator->(); }

      operator bool() const { return Base_ && Base_.index() < LastIndex_; }

      base_type const& base() const { return Base_; }
      base_type& base() { return Base_; }
      difference_type stride() const { return Stride_; }

   private:
      bool off_stride() const 
	 { return Base_ && Base_.index() < LastIndex_ && (Base_.index() - Start_) % Stride_ != 0; }

      base_type Base_;
      difference_type Start_;
      size_type Size_;
      difference_type Stride_;
      size_type LastIndex_;
}

// dense iterator, variable stride.

template <typename BaseIter>
class StrideIterator<BaseIter, tagVariable, vector_iterator_dense>
{
   public:
      typedef BaseIter base_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer   pointer;
      typedef typename base_type::category  category;

      StrideIterator(base_type const& Base, size_type Start, size_type Size, difference_type Stride,
		     difference_type Index = 0)
	 : Base_(Base), Size_(Size), Stride_(Stride), Index_(Index) { Base_ += Start; }

      StrideIterator& operator++()
      { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); ++Index_; return *this; }

      StrideIterator operator++(int) 
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); 
	   return StrideIterator(Base_, 0, Size_, Stride_, Index_++); }

      StrideIterator& operator--() 
         { DEBUG_CHECK(Index_ > 0)(Index_); --Index_; return *this; }

      StrideIterator operator--(int) 
         { DEBUG_CHECK(Index_ > 0)(Index_); 
	   return StrideIterator(Base_, 0, Size_, Stride_, Index_--); }

      StrideIterator& operator+=(difference_type n)
	 { DEBUG_CHECK(Index_ + n >= 0 && 
		       Index_ + n <= Size_)(Index_)(Size_);
	 Index_ += n; return *this; }

      StrideIterator& operator-=(difference_type n)
	 { DEBUG_CHECK(Index_ - n >= 0 && 
		       Index_ - n <= Size_)(Index_)(Size_);
	 Index_ -= n; return *this; }

      size_type index() const { return Index_; }

      reference operator*() const { return Base_[Index_ * Stride_]; }

      pointer operator->() const { return Base_ + Index_ * Stride_; }

      reference operator[](difference_type n) const 
      { return Base_[difference_type(Index_)+n]; }

      operator bool() const { return Index_ != Size_; }
	 
      base_type const& base() const { return Base_; }
      base_type& base() { return Base_; }
      difference_type stride() const { return Stride_; }

   private:
      base_type Base_;
      size_type Size_;
      difference_type Stride_;
      difference_type Index_;
};

// dense iterator, fixed stride.

template <typename BaseIter, int Stride_>
class StrideIterator<BaseIter, boost::mpl::int_<Stride_>, vector_iterator_dense>
{
   public:
      typedef BaseIter base_type;
      typedef typename base_type::reference reference;
      typedef typename base_type::pointer   pointer;
      typedef typename base_type::category  category;

      StrideIterator(base_type const& Base, size_type Start, size_type Size, 
		     difference_type Stride = Stride_,
		     difference_type Index = 0)
	 : Base_(Base), Size_(Size), Index_(Index) 
      { DEBUG_PRECONDITION_EQUAL(Stride, Stride_); Base_ += Start; }

      StrideIterator& operator++()
      { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); ++Index_; return *this; }

      StrideIterator operator++(int) 
         { DEBUG_CHECK(Index_ < Size_)(Index_)(Size_); 
	   return StrideIterator(Base_, 0, Size_, Stride_, Index_++); }

      StrideIterator& operator--() 
         { DEBUG_CHECK(Index_ > 0)(Index_); --Index_; return *this; }

      StrideIterator operator--(int) 
         { DEBUG_CHECK(Index_ > 0)(Index_); 
	   return StrideIterator(Base_, 0, Size_, Stride_, Index_--); }

      StrideIterator& operator+=(difference_type n)
	 { DEBUG_CHECK(Index_ + n >= 0 && 
		       Index_ + n <= Size_)(Index_)(Size_);
	 Index_ += n; return *this; }

      StrideIterator& operator-=(difference_type n)
	 { DEBUG_CHECK(Index_ - n >= 0 && 
		       Index_ - n <= Size_)(Index_)(Size_);
	 Index_ -= n; return *this; }

      size_type index() const { return Index_; }

      reference operator*() const { return Base_[Index_ * Stride_]; }

      pointer operator->() const { return Base_ + Index_ * Stride_; }

      reference operator[](difference_type n) const 
      { return Base_[difference_type(Index_)+n]; }

      operator bool() const { return Index_ != Size_; }
	 
      base_type const& base() const { return Base_; }
      base_type& base() { return Base_; }
      difference_type stride() const { return Stride_; }

      static difference_type const static_stride = Stride_;

   private:
      base_type Base_;
      size_type Size_;
      difference_type Index_;
};

} // namespace LinearAlgebra

#endif
