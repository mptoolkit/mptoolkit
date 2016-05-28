// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixiterators.h
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
  matrixiterators.h

  Created 2005-01-17 Ian McCulloch

  The beginnings of some adaptors to create coordinate iterators
  from matrix iterators.
*/

#if !defined(MATRIXITERATORS_H_H43Y598Y98WHHF8REHGLO83HLO)
#define MATRIXITERATORS_H_H43Y598Y98WHHF8REHGLO83HLO

#include "matrixinterface.h"

namespace LinearAlgebra
{

// IterAdaptorBase - implements everything except index1(), index2()
template <typename Iter, typename Derived>
class IterAdaptorBase
{
   private:
      BOOST_MPL_ASSERT_NOT((is_proxy_reference<Iter>));

     typedef Derived self_type;
     self_type const& self() const { return static_cast<self_type&>(*this); }
     self_type& self() { return static_cast<self_type&>(*this); }

   public:
     typedef Iter base_iterator;
     typedef typename base_iterator::value_type value_type;
     typedef typename base_iterator::reference  reference;
     typedef typename base_iterator::pointer    pointer;
     typedef typename base_iterator::category category;

      IterAdaptorBase() {}

      explicit IterAdaptorBase(Iter const& I) : I_(I) {}

      Derived& operator++() { ++I_; return this->self(); }
      Derived operator++(int);

      IterAdaptorBase& operator--() { --I_; return this->self(); }
      IterAdaptorBase operator--(int);

      IterAdaptorBase& operator+=(difference_type n)
	 { I_ += n; return this->self(); }

      IterAdaptorBase& operator-=(difference_type n)
	 { I_ -= n; return this->self(); }

      reference operator*() const { return *I_; }

      pointer operator->() const { return I_.operator->(); }

      reference operator[](difference_type n) const
	 { return I_[n]; }

      operator bool() const { return I_; }

   protected:
      size_type index() const { return I_.index(); }
      base_iterator const& base() const { return I_; }
      base_iterator& base() { return I_; }

   private:
      base_iterator I_;
};

template <typename Iter>
class MatrixIterFix1 : public IterAdaptorBase<Iter, MatrixIterFix1<Iter> >
{
   public:
      typedef IterAdaptorBase<Iter, MatrixIterFix1<Iter> > base_type;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::reference  reference;
      typedef typename base_type::pointer    pointer;

      MatrixIterFix1() {}

      MatrixIterFix1(Iter const& i, size_type Index1) : base_type(i), Index1_(Index1) {}

      size_type index1() const { return Index1_; }
      size_type index2() const { return this->index(); }

   private:
      size_type Index1_;
};

template <typename Iter>
class MatrixIterFix2 : public IterAdaptorBase<Iter, MatrixIterFix2<Iter> >
{
   public:
      typedef IterAdaptorBase<Iter, MatrixIterFix2<Iter> > base_type;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::reference  reference;
      typedef typename base_type::pointer    pointer;

      MatrixIterFix2() {}

      MatrixIterFix2(Iter const& i, size_type Index2) : base_type(i), Index2_(Index2) {}

      size_type index1() const { return this->index(); }
      size_type index2() const { return Index2_; }

   private:
      size_type Index2_;
};

//
// MatrixIterDiagonal
//

template <typename Iter>
class MatrixIterDiagonal : public IterAdaptorBase<Iter, MatrixIterDiagonal<Iter> >
{
   public:
      typedef IterAdaptorBase<Iter, MatrixIterDiagonal<Iter> > base_type;
      typedef typename base_type::value_type value_type;
      typedef typename base_type::reference  reference;
      typedef typename base_type::pointer    pointer;

      MatrixIterDiagonal() {}

      explicit MatrixIterDiagonal(Iter const& i) : base_type(i) {}

      size_type index1() const { return this->index(); }
      size_type index2() const { return this->index(); }
};

// adaptors

// MatrixDummyOuterIterator
// A trivial outer iterator for a coordinate type

template <typename Iter>
class MatrixDummyOuterIterator;

template <typename Iter>
struct Iterate<MatrixDummyOuterIterator<Iter>&>
{
   typedef Iter result_type;
   typedef MatrixDummyOuterIterator<Iter>& argument_type;
   result_type operator()(argument_type i) const { return i.iterate(); }
};

template <typename Iter>
struct Iterate<MatrixDummyOuterIterator<Iter> >
{
   typedef Iter result_type;
   typedef MatrixDummyOuterIterator<Iter> const& argument_type;
   result_type operator()(argument_type i) const { return i.iterate(); }
};

// MatrixOuterIterator<Iter, Orientation>
// An outer iterator for a vector of vectors, with a given orientation.
// Iter is an iterator over a vector of vectors.

template <typename Iter, typename Orientation>
class MatrixOuterIterator;

template <typename Iter, typename Orientation>
struct Iterate<MatrixOuterIterator<Iter, Orientation>&>
{
   typedef typename MatrixOuterIterator<Iter, Orientation>::iterator result_type;
   typedef MatrixOuterIterator<Iter, Orientation>& argument_type;
   result_type operator()(argument_type i) const { return i.iterate(); }
};

template <typename Iter, typename Orientation>
struct Iterate<MatrixOuterIterator<Iter, Orientation> >
{
   typedef typename MatrixOuterIterator<Iter, Orientation>::iterator result_type;
   typedef MatrixOuterIterator<Iter, Orientation> const& argument_type;
   result_type operator()(argument_type i) const { return i.iterate(); }
};

//
// MatrixDummyOuterIterator
//

template <typename Iter>
class MatrixDummyOuterIterator
{
   public:
      BOOST_MPL_ASSERT_NOT((is_proxy_reference<Iter>));

      typedef Iter iterator;

      typedef void reference;

      MatrixDummyOuterIterator() {}

      MatrixDummyOuterIterator(iterator const& i) : i_(i), end_(false) {}

      MatrixDummyOuterIterator(iterator const& i, bool end) : i_(i), end_(end) {}

      MatrixDummyOuterIterator& operator++() 
      { DEBUG_CHECK(!end_); end_ = true; return *this; }

      MatrixDummyOuterIterator operator++(int) 
      { DEBUG_CHECK(!end_); return MatrixDummyOuterIterator(i_, true); }

      operator bool() const { return !end_; }

      iterator iterate() const { return i_; }

   private:
      Iter i_;
      bool end_;
};

//
// MatrixOuterIterator
//

template <typename Iter, typename Orientation, typename OuterCategory>
class MatrixInnerIterator;

template <typename Iter>
class MatrixOuterIterator<Iter, RowMajor> : public Iter
{
   public:
      BOOST_MPL_ASSERT_NOT((is_proxy_reference<Iter>));

      typedef typename Iter::reference reference;
      typedef typename Iter::category category;
      typedef typename Iterate<typename reference_to_arg<reference>::type>::result_type 
      inner_vector_iterator;
      typedef MatrixInnerIterator<inner_vector_iterator, RowMajor, category> iterator;

      MatrixOuterIterator() {}

      MatrixOuterIterator(Iter const& i) : Iter(i) {}

      template <typename U>
      MatrixOuterIterator(U const& i) : Iter(i) {}

      size_type index1() const { return this->index(); }

      iterator iterate() const 
      { using LinearAlgebra::iterate; 
      return iterator(iterate(this->operator*()), this->index()); }
};

template <typename Iter>
class MatrixOuterIterator<Iter, ColMajor> : public Iter
{
   public:
      BOOST_MPL_ASSERT_NOT((is_proxy_reference<Iter>));

      typedef typename Iter::reference reference;
      typedef typename Iter::category category;
      typedef typename Iterate<typename reference_to_arg<reference>::type>::result_type 
      inner_vector_iterator;
      typedef MatrixInnerIterator<inner_vector_iterator, ColMajor, category> iterator;

      MatrixOuterIterator() {}

      MatrixOuterIterator(Iter const& i) : Iter(i) {}

      template <typename U>
      MatrixOuterIterator(U const& i) : Iter(i) {}

      size_type index2() const { return this->index(); }

      iterator iterate() const 
      { using LinearAlgebra::iterate; 
        return iterator(iterate(this->operator*()), this->index()); }
};

//
// MatrixInnerIterator
//

template <typename Iter, typename OuterCategory>
class MatrixInnerIterator<Iter, RowMajor, OuterCategory> 
   : public IterAdaptorBase<Iter,  MatrixInnerIterator<Iter, RowMajor, OuterCategory> >
{
   public:
      typedef IterAdaptorBase<Iter, 
			      MatrixInnerIterator<Iter, RowMajor, OuterCategory> > base_type;
      typedef typename get_matrix_category<OuterCategory, typename Iter::category>::type
          category;

      MatrixInnerIterator() {}

      MatrixInnerIterator(Iter const& i, size_type Index1) 
	 : base_type(i), Index1_(Index1) {}

      template <typename U>
      MatrixInnerIterator(U const& i, size_type Index1) : base_type(i), Index1_(Index1) {}

      size_type index1() const { return Index1_; }
      size_type index2() const { return this->index(); }
   
   private:
      size_type Index1_;
};

template <typename Iter, typename OuterCategory>
class MatrixInnerIterator<Iter, ColMajor, OuterCategory> 
   : public IterAdaptorBase<Iter, MatrixInnerIterator<Iter, ColMajor, OuterCategory> >
{
   public:
      BOOST_MPL_ASSERT_NOT((is_proxy_reference<Iter>));

      typedef IterAdaptorBase<Iter, 
			      MatrixInnerIterator<Iter, ColMajor, OuterCategory> > base_type;
      typedef typename get_matrix_category<OuterCategory, typename Iter::category>::type
          category;

      MatrixInnerIterator() {}

      MatrixInnerIterator(Iter const& i, size_type Index2) 
	 : base_type(i), Index2_(Index2) {}

      template <typename U>
      MatrixInnerIterator(U const& i, size_type Index2) : base_type(i), Index2_(Index2) {}

      size_type index1() const { return this->index(); }
      size_type index2() const { return Index2_; }
   
   private:
      size_type Index2_;
};

} // namespace LinearAlgebra

#endif
