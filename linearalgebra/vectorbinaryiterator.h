// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectorbinaryiterator.h
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
  vectorbinaryiterator.h

  Created 2005-01-10, Ian McCulloch
*/

#if !defined(VECTORBINARYITERATOR_H_CHHUIYTH879YO)
#define VECTORBINARYITERATOR_H_CHHUIYTH879YO

#include "vectortransformiterator.h"  // for operator_arrow_proxy

namespace LinearAlgebra
{


template <typename Iter1, typename Iter2, typename Func,
	  typename Cat1 = typename Iter1::category,
	  typename Cat2 = typename Iter2::category>
class VectorBinaryIterator;

template <typename Iter1, typename Iter2, typename Func>
class VectorBinaryIterator<Iter1, Iter2, Func, vector_iterator_dense, vector_iterator_dense>
{
   public:
      typedef Func functor_type;
      typedef typename make_value<typename functor_type::result_type>::type value_type;
      typedef typename functor_type::result_type reference;
      typedef operator_arrow_proxy<value_type> pointer;
      typedef vector_iterator_dense category;

      typedef Iter1 iterator1_type;
      typedef Iter2 iterator2_type;

      VectorBinaryIterator() {}

      explicit VectorBinaryIterator(iterator1_type const& i1, iterator2_type const& i2,
				    functor_type f = functor_type())
	 : i1_(i1), i2_(i2), f_(f) {}

      VectorBinaryIterator& operator++() 
	 { ++i1_; ++i2_; return *this; }

      VectorBinaryIterator& operator++(int) 
	 { return VectorBinaryIterator(i1_++, i2_++, f_); }

      VectorBinaryIterator& operator--() 
	 { --i1_; --i2_; return *this; }

      VectorBinaryIterator& operator--(int) 
	 { return VectorBinaryIterator(i1_--, i2_--, f_); }

      VectorBinaryIterator& operator+=(difference_type n)
	 { i1_ += n; i2_ += n; return *this; }

      VectorBinaryIterator& operator-=(difference_type n)
	 { i1_ -= n; i2_ -= n; return *this; }

      size_type index() const 
      { DEBUG_PRECONDITION(bool(i1_));
        DEBUG_CHECK_EQUAL(i1_.index(), i2_.index()); 
        return i1_.index(); }

      reference operator*() const { return f_(*i1_, *i2_); }

      pointer operator->() const { return pointer(&f_(*i1_, *i2_)); }

      operator bool() const { DEBUG_CHECK_EQUAL(bool(i1_), bool(i2_)); return i1_; }

      iterator1_type const& iterator1() const { return i1_; }
      iterator2_type const& iterator2() const { return i2_; }

      functor_type const& functor() const { return f_; }

      // derived concepts
      size_type size() const { return i1_.size(); }
      reference operator[](difference_type n) const { return f_(i1_[n], i2_[n]); }

   private:
      iterator1_type i1_;
      iterator2_type i2_;
      functor_type f_;
};
   
} // namespace LinearAlgebra

#endif
