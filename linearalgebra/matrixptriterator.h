// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixptriterator.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
/* -*- C++ -*- $Id$

  matrixptriterator.h

*/

#if !defined(MATRIXPTRITERATOR_H_CJKY89YU89P5Y89)
#define MATRIXPTRITERATOR_H_CJKY89YU89P5Y89

#include "vecptriterator.h"
#include "vectormemproxy.h"
#include "vector.h"
#include "matrixinterface.h"
#include "matrixiterators.h"
#include "matrixoperationsbase.h"

namespace LinearAlgebra
{

// TODO: eventually, this should be templated on the size & stride

// TODO: fix up orientation
template <typename Scalar, typename Orientation>
class MatrixPtrIterator;

template <typename Scalar>
class MatrixPtrIterator<Scalar, RowMajor>
{
   public:
      BOOST_MPL_ASSERT_NOT((boost::is_reference<Scalar>));

      typedef VectorMemProxy<Scalar, tagVariable> reference;
      typedef typename make_value<reference>::type value_type;

      typedef VecPtrIterator<Scalar, tagVariable> basic_iterator;
      typedef MatrixIterFix1<basic_iterator> iterator;
      typedef VecPtrIterator<Scalar const, tagVariable> const_basic_iterator;
      typedef MatrixIterFix1<const_basic_iterator> const_iterator;
      typedef operator_arrow_proxy<reference> pointer;

      typedef vector_iterator_dense category;

      MatrixPtrIterator(Scalar* Base,
                        size_type Size1, difference_type Stride1,
                        size_type Size2, difference_type Stride2,
                        size_type OuterIndex = 0)
         : Base_(Base), Size1_(Size1), Stride1_(Stride1),
         Size2_(Size2), Stride2_(Stride2),
         OuterIndex_(OuterIndex)
         {}

      template <typename OtherScalar>
      MatrixPtrIterator(MatrixPtrIterator<OtherScalar, RowMajor> const& Other)
         : Base_(Other.base()), Size1_(Other.size1()), Stride1_(Other.stride1()),
         Size2_(Other.size2()), Stride2_(Other.stride2()), OuterIndex_(Other.index())
        {}
      //{ TRACE(Base_); }

      MatrixPtrIterator& operator++()
         { DEBUG_CHECK(OuterIndex_ < Size1_)(OuterIndex_)(Size1_);
         ++OuterIndex_; return *this; }

      MatrixPtrIterator operator++(int)
         { DEBUG_CHECK(OuterIndex_ < Size1_)(OuterIndex_)(Size1_);
           return MatrixPtrIterator(Base_, Size1_, Stride1_,
                                    Size2_, Stride2_, OuterIndex_++); }

      MatrixPtrIterator& operator--()
         { DEBUG_CHECK(OuterIndex_ > 0)(OuterIndex_);
         --OuterIndex_; return *this; }

      MatrixPtrIterator operator--(int)
         { DEBUG_CHECK(OuterIndex_ > 0)(OuterIndex_);
           return MatrixPtrIterator(Base_, Size1_, Stride1_,
                                    Size2_, Stride2_, OuterIndex_--); }

      MatrixPtrIterator& operator+=(difference_type n)
         { DEBUG_CHECK(difference_type(OuterIndex_) + n >= 0 &&
                       difference_type(OuterIndex_) + n <= Size1_)
                       (OuterIndex_)(Size1_);
         OuterIndex_ += n; return *this; }

      MatrixPtrIterator& operator-=(difference_type n)
         { DEBUG_CHECK(difference_type(OuterIndex_) - n >= 0 &&
                       difference_type(OuterIndex_) - n <= Size1_)
                       (OuterIndex_)(Size1_);
         OuterIndex_ -= n; return *this; }

      size_type index() const { return OuterIndex_; }

      reference operator*() const
         { return reference(this->ptr(), Size2_, Stride2_); }

      reference operator[](difference_type n) const
         { return reference(this->ptr(n), Size2_, Stride2_); }

      operator bool() const { return OuterIndex_ != Size1_; }

      Scalar* ptr() const { return Base_ + difference_type(OuterIndex_) * Stride1_; }
      Scalar* ptr(difference_type n) const
         { return Base_ + (difference_type(OuterIndex_) + n) * Stride1_; }

      size_type size() const { return Size1_ - OuterIndex_; }

      size_type size1() const { return Size1_; }
      difference_type stride1() const { return Stride1_; }

      size_type size2() const { return Size2_; }
      difference_type stride2() const { return Stride2_; }

      difference_type outer_stride() const { return Stride1_; }
      size_type inner_size() const { return Size2_; }
      difference_type inner_stride() const { return Stride2_; }

      iterator iterate() const
         { return iterator(basic_iterator(this->ptr(), Size2_, 0, Stride2_),
                           OuterIndex_); }

      const_iterator const_iterate() const
         { return const_iterator(const_basic_iterator(this->ptr(), Size2_, 0, Stride2_),
                                 OuterIndex_); }

      Scalar* base() const { return Base_; }

   private:
      Scalar* Base_;
      size_type Size1_;
      difference_type Stride1_;
      size_type Size2_;
      difference_type Stride2_;
      size_type OuterIndex_;
};

template <typename Scalar>
class MatrixPtrIterator<Scalar, ColMajor>
{
   public:
      typedef VectorMemProxy<Scalar, tagVariable> reference;
      typedef typename make_value<reference>::type value_type;

      typedef VecPtrIterator<Scalar, tagVariable> basic_iterator;
      typedef MatrixIterFix2<basic_iterator> iterator;
      typedef VecPtrIterator<Scalar const, tagVariable> const_basic_iterator;
      typedef MatrixIterFix2<const_basic_iterator> const_iterator;
      typedef operator_arrow_proxy<reference> pointer;

      typedef vector_iterator_dense category;

      MatrixPtrIterator(Scalar* Base,
                        size_type Size1, difference_type Stride1,
                        size_type Size2, difference_type Stride2,
                        size_type OuterIndex = 0)
         : Base_(Base), Size1_(Size1), Stride1_(Stride1),
         Size2_(Size2), Stride2_(Stride2),
         OuterIndex_(OuterIndex)
         {}

      template <typename OtherScalar>
      MatrixPtrIterator(MatrixPtrIterator<OtherScalar, ColMajor> const& Other)
         : Base_(Other.base()), Size1_(Other.size1()), Stride1_(Other.stride1()),
         Size2_(Other.size2()), Stride2_(Other.stride2()), OuterIndex_(Other.index()) {}

      MatrixPtrIterator& operator++()
         { DEBUG_CHECK(OuterIndex_ < Size2_)(OuterIndex_)(Size2_);
         ++OuterIndex_; return *this; }

      MatrixPtrIterator operator++(int)
         { DEBUG_CHECK(OuterIndex_ < Size2_)(OuterIndex_)(Size2_);
           return MatrixPtrIterator(Base_, Size1_, Stride1_,
                                    Size2_, Stride2_, OuterIndex_++); }

      MatrixPtrIterator& operator--()
         { DEBUG_CHECK(OuterIndex_ > 0)(OuterIndex_);
         --OuterIndex_; return *this; }

      MatrixPtrIterator operator--(int)
         { DEBUG_CHECK(OuterIndex_ > 0)(OuterIndex_);
           return MatrixPtrIterator(Base_, Size1_, Stride1_,
                                    Size2_, Stride2_, OuterIndex_--); }

      MatrixPtrIterator& operator+=(difference_type n)
         { DEBUG_CHECK(difference_type(OuterIndex_) + n >= 0 &&
                       difference_type(OuterIndex_) + n <= Size2_)
                       (OuterIndex_)(Size2_);
         OuterIndex_ += n; return *this; }

      MatrixPtrIterator& operator-=(difference_type n)
         { DEBUG_CHECK(difference_type(OuterIndex_) - n >= 0 &&
                       difference_type(OuterIndex_) - n <= Size2_)
                       (OuterIndex_)(Size2_);
         OuterIndex_ -= n; return *this; }

      size_type index() const { return OuterIndex_; }

      reference operator*() const
         { return reference(this->ptr(), Size1_, Stride1_); }

      reference operator[](difference_type n) const
         { return reference(this->ptr(n), Size1_, Stride1_); }

      operator bool() const { return OuterIndex_ != Size2_; }

      Scalar* ptr() const { return Base_ + difference_type(OuterIndex_) * Stride2_; }
      Scalar* ptr(difference_type n) const
         { return Base_ + (difference_type(OuterIndex_) + n) * Stride2_; }

      size_type size() const { return Size2_ - OuterIndex_; }

      size_type size1() const { return Size1_; }
      difference_type stride1() const { return Stride1_; }

      size_type size2() const { return Size2_; }
      difference_type stride2() const { return Stride2_; }

      difference_type outer_stride() const { return Stride2_; }
      size_type inner_size() const { return Size1_; }
      difference_type inner_stride() const { return Stride1_; }

      iterator iterate() const
         { return iterator(basic_iterator(this->ptr(), Size1_, 0, Stride1_),
                           OuterIndex_); }

      const_iterator const_iterate() const
         { return const_iterator(const_basic_iterator(this->ptr(), Size1_, 0, Stride1_),
                                 OuterIndex_); }

      Scalar* base() const { return Base_; }

   private:
      Scalar* Base_;
      size_type Size1_;
      difference_type Stride1_;
      size_type Size2_;
      difference_type Stride2_;
      size_type OuterIndex_;
};

// iterators

template <typename Scalar, typename Orient>
struct Iterate<MatrixPtrIterator<Scalar, Orient>&>
{
   typedef typename MatrixPtrIterator<Scalar, Orient>::iterator result_type;
   typedef MatrixPtrIterator<Scalar, Orient> argument_type;
   result_type operator()(argument_type const& x) const
   { return x.iterate(); }
};

template <typename Scalar, typename Orient>
struct Iterate<MatrixPtrIterator<Scalar, Orient> >
{
   typedef typename MatrixPtrIterator<Scalar, Orient>::const_iterator result_type;
   typedef MatrixPtrIterator<Scalar, Orient> argument_type;
   result_type operator()(argument_type const& x) const
   { return x.const_iterate(); }
};


} // namespace LinearAlgebra

#endif
