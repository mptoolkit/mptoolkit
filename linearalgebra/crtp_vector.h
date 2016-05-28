// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/crtp_vector.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  crtp_vector.h

  Curiously-recursive base classes to define some useful member functions.

  Created 2005-02-28 Ian McCulloch
*/

#if !defined(CRTP_VECTOR_H_HDSUIHUY789YT78YWE9)
#define CRTP_VECTOR_H_HDSUIHUY789YT78YWE9

#include "vectoroperationsbase.h"
#include <boost/type_traits.hpp>

namespace LinearAlgebra
{

template <typename Derived>
class VectorBase
{
   private:
      typedef Derived derived_type;

      derived_type& as_derived() { return static_cast<derived_type&>(*this); }

      derived_type const& as_derived() const
      { return static_cast<derived_type const&>(*this); }

   public:

      //      typedef typename interface<derived_type>::value_type value_type;

      // size

      typename Size<derived_type>::result_type
      size() const { return Size<derived_type>()(this->as_derived()); }

      // resize

      template <typename Int>
      typename boost::enable_if<boost::is_convertible<Int, size_type>,
         Resize<derived_type&> >::type::result_type
      resize(Int n)
      {
         return Resize<derived_type&>()(this->as_derived(), n);
      }

      // operator +=

      template <typename RHS>
      typename boost::enable_if<is_defined<AddCopy<derived_type&, RHS> >, VectorBase&>::type
      operator+=(RHS const& x)
      {
         add_copy(this->as_derived(), x);
         return *this;
      }

#if 0
      template <typename RHS>
      typename AddCopy<derived_type&, value_with_zero<RHS> >::result_type
      operator+=(value_with_zero<RHS> const& x)
      {
         add_copy(this->as_derived(), x);
         return *this;
      }
#endif

      // operator -=

      template <typename RHS>
      typename boost::enable_if<is_defined<SubtractCopy<derived_type&, RHS> >, VectorBase&>::type
      operator-=(RHS const& x)
      {
         subtract_copy(this->as_derived(), x);
         return *this;
      }

#if 0
      template <typename RHS>
      typename SubtractCopy<derived_type&, value_with_zero<RHS> >::result_type
      operator-=(value_with_zero<RHS> const& x)
      {
         return subtract(this->as_derived(), x);
      }
#endif

      // operator *=

      template <typename RHS>
      typename Multiply<derived_type&, RHS>::result_type
      operator*=(RHS const& x)
      {
         return multiply(this->as_derived(), x);
      }

      // operator []

      template <typename RHS>
      typename boost::disable_if<is_const_proxy<derived_type>, 
                                 typename VectorBracket<derived_type&, RHS>::result_type>::type
      operator[](RHS const& x)
      {
         return VectorBracket<derived_type&, RHS>()(this->as_derived(), x);
      }

      // operator [] const

      template <typename RHS>
      typename VectorBracket<derived_type, RHS>::result_type
      operator[](RHS const& x) const
      {
         return VectorBracket<derived_type, RHS>()(this->as_derived(), x);
      }
};

} // namespace LinearAlgebra

#endif
