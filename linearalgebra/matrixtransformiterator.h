// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/matrixtransformiterator.h
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
  matrixtransformiterator.h

  Created 2005-01-24 Ian McCulloch
*/

#if !defined(MPTOOLKIT_LINEARALGEBRA_MATRIXTRANSFORMITERATOR_H)
#define MPTOOLKIT_LINEARALGEBRA_MATRIXTRANSFORMITERATOR_H

#include "vectortransformiterator.h"
#include "matrixiterators.h"

namespace LinearAlgebra
{

// TODO: some kind of iterator facade is needed here

template <typename Base, typename Func>
class MatrixTransformInnerIterator
{
   public:
      typedef typename Func::result_type reference;
      typedef typename make_value<reference>::type value_type;
      typedef typename boost::mpl::if_<
	 boost::is_reference<reference>,
	 typename boost::remove_reference<reference>::type*,
	 operator_arrow_proxy<value_type> >::type pointer;
      typedef typename Base::category category;

      MatrixTransformInnerIterator() {}

      explicit MatrixTransformInnerIterator(Base const& base)
	 : base_(base), f_() {}

      MatrixTransformInnerIterator(Base const& base, Func const& f)
	 : base_(base), f_(f) {}

      MatrixTransformInnerIterator& operator++() 
	 { ++base_; return *this; }

      MatrixTransformInnerIterator& operator++(int) 
	 { return MatrixTransformInnerIterator(base_++, f_); }

      MatrixTransformInnerIterator& operator--() 
	 { --base_; return *this; }

      MatrixTransformInnerIterator& operator--(int) 
	 { return MatrixTransformInnerIterator(base_--, f_); }

      MatrixTransformInnerIterator& operator+=(difference_type n)
	 { base_ += n; return *this; }

      MatrixTransformInnerIterator& operator-=(difference_type n)
	 { base_ -= n; return *this; }

      size_type index1() const { return base_.index1(); }
      size_type index2() const { return base_.index2(); }

      reference operator*() const { return f_(*base_); }

      pointer operator->() const { return pointer(&f_(*base_)); }

      operator bool() const { return base_; }

      Base& base() { return base_; }
      Base const& base() const { return base_; }

      // derived concepts
   //      size_type size() const { return base_.size(); }
      reference operator[](difference_type n) const { return f_(base_[n]); }
      reference operator()(size_type n) const { return f_(base_(n)); }
      bool element_exists(size_type n) const { return base_.element_exists(n); }

   private:
      Base base_;
      Func f_;
};

template <typename Base, typename Func>
class MatrixTransformOuterIterator
{
   public:
      // if F allows a mutable parameter, then we want to use the
      // mutable version of Transform (ie. Transform<T&, Func>).
      typedef typename boost::mpl::if_<
	 is_mutable_proxy_reference<typename Func::argument_type>,
	 typename boost::add_reference<typename Base::reference>::type,
	 typename basic_type<typename Base::reference>::type
      >::type base_arg_type;

      typedef typename Transform<base_arg_type, Func>::result_type reference;

   //   typedef typename boost::mpl::print<base_arg_type>::type d;
   //   typedef typename boost::mpl::print<reference>::type d2;

      //typedef typename Transform<typename Base::reference, Func>::result_type reference;
      typedef typename make_value<reference>::type value_type;
      typedef typename boost::mpl::if_<
	 boost::is_reference<reference>,
	 typename boost::remove_reference<reference>::type*,
	 operator_arrow_proxy<value_type> >::type pointer;
      typedef typename Base::category category;

      typedef typename iterator<Base>::type base_iterator;
      typedef MatrixTransformInnerIterator<base_iterator, Func> iterator;

      MatrixTransformOuterIterator() {}

      explicit MatrixTransformOuterIterator(Base const& base)
	 : base_(base), f_() {}

      MatrixTransformOuterIterator(Base const& base, Func const& f)
	 : base_(base), f_(f) {}

      MatrixTransformOuterIterator& operator++() 
	 { ++base_; return *this; }

      MatrixTransformOuterIterator& operator++(int) 
	 { return MatrixTransformOuterIterator(base_++, f_); }

      MatrixTransformOuterIterator& operator--() 
	 { --base_; return *this; }

      MatrixTransformOuterIterator& operator--(int) 
	 { return MatrixTransformOuterIterator(base_--, f_); }

      MatrixTransformOuterIterator& operator+=(difference_type n)
	 { base_ += n; return *this; }

      MatrixTransformOuterIterator& operator-=(difference_type n)
	 { base_ -= n; return *this; }

      size_type index() const { return base_.index(); }

      reference operator*() const { return transform(*base_, f_); }

   //      pointer operator->() const { return pointer(&transform(*base_, f_)); }

      operator bool() const { return base_; }

      Base& base() { return base_; }
      Base const& base() const { return base_; }

      iterator iterate() const 
      { using LinearAlgebra::iterate; return iterator(iterate(base_), f_); }

      // derived concepts
      // size_type size() const { return base_.size(); }
      reference operator[](difference_type n) const { return transform(base_[n], f_); }
      reference operator()(size_type n) const { return transform(base_(n), f_); }
      bool element_exists(size_type n) const { return base_.element_exists(n); }

      Func const& func() const { return f_; }

   private:
      Base base_;
      Func f_;
};

// iterators

template <typename Base, typename Func>
struct Iterate<MatrixTransformOuterIterator<Base, Func>&>
{
   typedef typename MatrixTransformOuterIterator<Base, Func>::iterator result_type;
   typedef MatrixTransformOuterIterator<Base, Func>& argument_type;
   result_type operator()(argument_type x) const 
   { return result_type(iterate(x.base()), x.func()); }
};

template <typename Base, typename Func>
struct Iterate<MatrixTransformOuterIterator<Base, Func> >
{
   typedef typename MatrixTransformOuterIterator<Base, Func>::iterator result_type;
   typedef MatrixTransformOuterIterator<Base, Func> const& argument_type;
   result_type operator()(argument_type x) const 
   { return result_type(iterate(x.base()), x.func()); }
};

} // namespace LinearAlgebra

#endif
