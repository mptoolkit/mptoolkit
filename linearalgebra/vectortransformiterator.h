// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectortransformiterator.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
/* -*- C++ -*- $Id$

  vectortransformiterator.h

  abstracted out of vectortransform.h, 2005-01-09, Ian McCulloch
*/

#if !defined(VECTORTRANSFORMITERATOR_H_CHHUIYTH879YO)
#define VECTORTRANSFORMITERATOR_H_CHHUIYTH879YO

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>

namespace LinearAlgebra
{

// operator_arrow_proxy, 'borrowed' from boost::iterator_facade,
// slightly modified
template <class T>
struct operator_arrow_proxy
{
   template <typename U>
   operator_arrow_proxy(U const* px) : m_value(*px) { TRACE("ctor"); }
   ~operator_arrow_proxy() { TRACE("dtor"); }
   operator_arrow_proxy(operator_arrow_proxy const& x) : m_value(x.m_value) { TRACE("copy"); }
   const T* operator->() const { return &m_value; }
   // This function is needed for MWCW and BCC, which won't call operator->
   // again automatically per 13.3.1.2 para 8
   operator const T*() const { return &m_value; }
   T m_value;
};

template <typename Base, typename Func>
class VectorTransformIterator
{
   public:
      typedef typename Func::result_type reference;
      typedef typename make_value<reference>::type value_type;
      typedef typename boost::mpl::if_<
	 boost::is_reference<reference>,
	 typename boost::remove_reference<reference>::type*,
	 operator_arrow_proxy<value_type> >::type pointer;
      typedef typename Base::category category;

      VectorTransformIterator() {}

      explicit VectorTransformIterator(Base const& base)
	 : base_(base), f_() {}

      VectorTransformIterator(Base const& base, Func const& f)
	 : base_(base), f_(f) {}

      VectorTransformIterator& operator++() 
	 { ++base_; return *this; }

      VectorTransformIterator& operator++(int) 
	 { return VectorTransformIterator(base_++, f_); }

      VectorTransformIterator& operator--() 
	 { --base_; return *this; }

      VectorTransformIterator& operator--(int) 
	 { return VectorTransformIterator(base_--, f_); }

      VectorTransformIterator& operator+=(difference_type n)
	 { base_ += n; return *this; }

      VectorTransformIterator& operator-=(difference_type n)
	 { base_ -= n; return *this; }

      size_type index() const { return base_.index(); }


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

} // namespace LinearAlgebra

#endif
