// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/crtp_matrix.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ian@qusim.net>
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
  crtp_matrix.h

  Curiously-recursive base classes to define some useful member functions.

  Created 2005-03-03 Ian McCulloch
*/

#if !defined(CRTP_MATRIX_H_SJDLKHJGURIHY8979843)
#define CRTP_MATRIX_H_SJDLKHJGURIHY8979843

#include "matrixoperationsbase.h"
#include <boost/type_traits.hpp>

namespace LinearAlgebra
{

template <typename Derived>
class MatrixBase
{
   private:
      typedef Derived derived_type;

      derived_type& as_derived() { return static_cast<derived_type&>(*this); }

      derived_type const& as_derived() const
      { return static_cast<derived_type const&>(*this); }

   public:

      //      typedef typename interface<derived_type>::value_type value_type;

      // size1

      typename Size1<derived_type>::result_type
      size1() const { return Size1<derived_type>()(this->as_derived()); }

      // size2

      typename Size2<derived_type>::result_type
      size2() const { return Size2<derived_type>()(this->as_derived()); }

      // resize

      template <typename Int1, typename Int2>
      typename boost::enable_if<
         boost::mpl::and_<
            boost::is_convertible<Int1, size_type>
          , boost::is_convertible<Int2, size_type>
         >
       , Resize<derived_type&>
      >::type::result_type
      resize(Int1 r, Int2 c)
      {
         return Resize<derived_type&>()(this->as_derived(), r, c);
      }

      // fill

      template <typename U>
      typename Fill<derived_type&, U>::result_type
      fill(U const& x)
      {
         return  Fill<derived_type&, U>(this->as_derived(), x);
      }

      // operator +=

      template <typename RHS>
      derived_type&
      //typename AddCopy<derived_type&, RHS>::result_type
      operator+=(RHS const& x)
      {
         Add<derived_type&, RHS>().operator()(this->as_derived(), x);
         return this->as_derived();
         //         return add_copy(this->as_derived(), x);
      }

      template <typename RHS>
      typename AddCopy<derived_type&, value_with_zero<RHS>>::result_type
      operator+=(value_with_zero<RHS> const& x)
      {
         if (!x.is_zero())
            add_copy(this->as_derived(), x.get());
         return *this;
      }

      // operator -=

      template <typename RHS>
      typename SubtractCopy<derived_type&, RHS>::result_type
      operator-=(RHS const& x)
      {
         return subtract_copy(this->as_derived(), x);
      }

      template <typename RHS>
      typename SubtractCopy<derived_type&, value_with_zero<RHS> >::result_type
      operator-=(value_with_zero<RHS> const& x)
      {
         return subtract_copy(this->as_derived(), x);
      }

      // operator *=

      template <typename RHS>
      typename Multiply<derived_type&, RHS>::result_type
      operator*=(RHS const& x)
      {
         return multiply(this->as_derived(), x);
      }

      // operator ()

      template <typename RHS1, typename RHS2>
      typename MatrixBracket<derived_type&, RHS1, RHS2>::result_type
      operator()(RHS1 const& x, RHS2 const& y)
      {
         return MatrixBracket<derived_type&, RHS1, RHS2>()(this->as_derived(), x, y);
      }

      // operator () const

      template <typename RHS1, typename RHS2>
      typename MatrixBracket<derived_type, RHS1, RHS2>::result_type
      operator()(RHS1 const& x, RHS2 const& y) const
      {
         //typedef typename boost::mpl::print<typename MatrixBracket<derived_type, RHS1, RHS2>::result_type>::type dummy;
         return MatrixBracket<derived_type, RHS1, RHS2>()(this->as_derived(), x, y);
      }
};

} // namespace LinearAlgebra

#endif
