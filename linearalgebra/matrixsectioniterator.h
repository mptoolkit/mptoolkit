// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/matrixsectioniterator.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
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

/*
  matrixsectioniterator.h

  Iterate over a section of a matrix, defined by row & column
  index vectors.

  Created 2005-02-28 Ian McCulloch
*/

#if !defined(MATRIXSECTIONITERATOR_H_JHRFY897Y789Y89YP)
#define MATRIXSECTIONITERATOR_H_JHRFY897Y789Y89YP

#include "matrixiterators.h"

namespace LinearAlgebra
{

template <typename Base, typename Orientation,
          typename OuterIndexIter, typename InnerIndexRef>
class MatrixSectionOuterIterator;

template <typename Base, typename RowIndexIter, typename ColIndexRef>
class MatrixSectionOuterIterator<Base, RowMajor, RowIndexIter, ColIndexRef>
{
   public:
      BOOST_MPL_ASSERT((is_const_proxy_reference<ColIndexRef>));

      typedef typename Base::value_type value_type;
      typedef typename Base::reference base_ref_type;
      typedef IndirectVector<base_ref_type, ColIndexRef> reference;
      typedef typename iterator<reference>::type base_iter_type;
      typedef MatrixIterFix1<base_iter_type> iterator;
      typedef operator_arrow_proxy<reference> pointer;
      typedef typename RowIndexIter::category category;

      MatrixSectionOuterIterator() {}

      MatrixSectionOuterIterator(Base b, RowIndexIter r, ColIndexRef c,
                                 size_type RowIndex = 0)
         : Base_(b), RowIndexIter_(r), ColIndex_(c), RowIndex_(RowIndex) {}

      template <typename OtherBase>
         MatrixSectionOuterIterator(MatrixSectionOuterIterator<
                                    OtherBase,
                                    RowMajor,
                                    RowIndexIter,
                                    ColIndexRef> const& Other)
         : Base_(Other.base()), RowIndexIter_(Other.row_index_iter()),
         ColIndex_(Other.col_index()) {}

      iterator iterate() const
   { using LinearAlgebra::iterate; return iterator(iterate(this->operator*()), this->index()); }

      MatrixSectionOuterIterator& operator++() { ++RowIndexIter_; ++RowIndex_; return *this; }

      MatrixSectionOuterIterator& operator++(int)
         { return MatrixSectionOuterIterator(Base_, RowIndexIter_++, ColIndex_, RowIndex_++); }

      MatrixSectionOuterIterator& operator+=(size_type n)
         { RowIndexIter_ += n; RowIndex_ += n; return *this; }

      reference operator*() const
         { return reference(Base_[*RowIndexIter_], ColIndex_); }

      reference operator[](difference_type n) const
         { return reference(Base_[RowIndexIter_[n]], ColIndex_); }

      pointer operator->() const { return pointer(&(this->operator*())); }

      operator bool() const { return RowIndexIter_; }

      size_type index() const { return RowIndexIter_.index(); }

      Base base() const { return Base_; }
      RowIndexIter row_index_iter() const { return RowIndexIter_; }
      ColIndexRef col_index() const { return ColIndex_; }

   private:
      Base Base_;
      RowIndexIter RowIndexIter_;
      ColIndexRef ColIndex_;
      size_type RowIndex_;
};

template <typename Base, typename ColIndexIter, typename RowIndexRef>
class MatrixSectionOuterIterator<Base, ColMajor, ColIndexIter, RowIndexRef>
{
   public:
      BOOST_MPL_ASSERT((is_const_proxy_reference<RowIndexRef>));

      typedef typename Base::value_type value_type;
      typedef typename Base::reference base_ref_type;
      typedef IndirectVector<base_ref_type, RowIndexRef> reference;
      typedef typename iterator<reference>::type base_iter_type;
      typedef MatrixIterFix2<base_iter_type> iterator;
      typedef operator_arrow_proxy<reference> pointer;
      typedef typename ColIndexIter::category category;

      MatrixSectionOuterIterator() {}

      MatrixSectionOuterIterator(Base b, ColIndexIter r, RowIndexRef c)
         : Base_(b), ColIndexIter_(r), ColIndex_(c) {}

      template <typename OtherBase>
         MatrixSectionOuterIterator(MatrixSectionOuterIterator<
                                    OtherBase,
                                    ColMajor,
                                    ColIndexIter,
                                    RowIndexRef> const& Other)
         : Base_(Other.base()), ColIndexIter_(Other.row_index_iter()),
         ColIndex_(Other.col_index()) {}

      iterator iterate() const
   { using LinearAlgebra::iterate; return iterator(iterate(this->operator*()), this->index()); }

      MatrixSectionOuterIterator& operator++() { ++ColIndexIter_; return *this; }

      MatrixSectionOuterIterator& operator++(int)
         { return MatrixSectionOuterIterator(Base_, ColIndexIter_++, ColIndex_); }

      MatrixSectionOuterIterator& operator+=(size_type n)
         { ColIndexIter_ += n; return *this; }

      reference operator*() const
         { return reference(Base_[*ColIndexIter_], ColIndex_); }

      reference operator[](difference_type n) const
         { return reference(Base_[ColIndexIter_[n]], ColIndex_); }

      pointer operator->() const { return pointer(&(this->operator*())); }

      operator bool() const { return ColIndexIter_; }

      size_type index() const { return ColIndexIter_.index(); }

      Base base() const { return Base_; }
      ColIndexIter row_index_iter() const { return ColIndexIter_; }
      RowIndexRef col_index() const { return ColIndex_; }

   private:
      Base Base_;
      ColIndexIter ColIndexIter_;
      RowIndexRef ColIndex_;
};

template <typename Base, typename Orient, typename OuterIndexIter, typename InnerIndexRef>
struct Iterate<MatrixSectionOuterIterator<Base, Orient, OuterIndexIter, InnerIndexRef> >
{
   typedef MatrixSectionOuterIterator<Base, Orient, OuterIndexIter, InnerIndexRef> base;
   typedef typename base::iterator result_type;
   typedef base const& argument_type;
   result_type operator()(argument_type i) const
   {
      return i.iterate();
   }
};

template <typename Base, typename Orient, typename OuterIndexIter, typename InnerIndexRef>
struct Iterate<MatrixSectionOuterIterator<Base, Orient, OuterIndexIter, InnerIndexRef>&>
{
   typedef MatrixSectionOuterIterator<Base, Orient, OuterIndexIter, InnerIndexRef> base;
   typedef typename base::iterator result_type;
   typedef base const& argument_type;
   result_type operator()(argument_type i) const
   {
      return i.iterate();
   }
};

} // namespace LinearAlgebra

#endif
