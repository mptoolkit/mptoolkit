// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/noalias.h
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
/* -*- C++ -*- $Id$

  noalias.h

  defines NoAliasProxy and noalias function, to declare that the left hand side
  of an expression does not alias with the right hand side.

  created 2005-01-07 Ian McCulloch
*/

#if !defined(NOALIAS_H_JKSFDHEJIRH4389UT89YP)
#define NOALIAS_H_JKSFDHEJIRH4389UT89YP

#include <boost/mpl/assert.hpp>

namespace LinearAlgebra
{

namespace Private
{

template <typename T, typename Enable = void>
struct is_declared_noalias : boost::mpl::false_ {};

template <typename T>
struct is_declared_noalias<T,
    typename boost::enable_if<exists<typename T::noalias> >::type>
: boost::mpl::true_ {};

} // namespace Private

template <typename T>
struct is_noalias : Private::is_declared_noalias<T> {};

template <typename T>
inline
T& strip_noalias(T& x)
{
   return x;
}

template <typename T>
inline
T const& strip_noalias(T const& x)
{
   return x;
}

template <typename RefType,
          typename RInterface = typename interface<typename make_value<RefType>::type>::type>
class NoAliasProxy
{
   public:
      BOOST_MPL_ASSERT((is_proxy_reference<RefType>));

      typedef boost::mpl::true_ proxy;

      typedef boost::mpl::true_ noalias;

      typedef typename boost::remove_const<
         typename boost::remove_reference<RefType>::type>::type value_type;

      typedef RefType reference;
      typedef typename make_const_reference<reference>::type const_reference;

      explicit NoAliasProxy(reference Value) : Value_(Value) {}

      // compiler-synthesized copy-ctor is OK

      // is this also the copy-assignment operator?
      template <typename U>
      NoAliasProxy& operator=(U const& x)
      {
         assign(Value_, x);
         return *this;
      }

      // just in case...
      NoAliasProxy& operator=(NoAliasProxy const& x)
      {
         assign(Value_, x.value());
         return *this;
      }

      // to allow assignment with temporary proxies,
      // computed assignment operators must be members here,
      // or see comment below for an alternative

      template <typename U>
      NoAliasProxy& operator+=(U const& x)
      {
         add(Value_, x);
         return *this;
      }

      template <typename U>
      NoAliasProxy& operator-=(U const& x)
      {
         subtract(Value_, x);
         return *this;
      }

      reference value() { return Value_; }
      const_reference value() const { return Value_; }

      operator NoAliasProxy<const_reference>() const
      { return NoAliasProxy<const_reference>(Value_); }

   private:
      reference Value_;
};

template <typename RefType, typename Rv, typename Ri>
class NoAliasProxy<RefType, ANY_VECTOR(Rv, Ri)>
{
   public:
      BOOST_MPL_ASSERT((is_proxy_reference<RefType>));

      typedef boost::mpl::true_ proxy;

      typedef typename boost::remove_const<typename boost::remove_reference<RefType>::type>::type value_type;

      typedef RefType reference;
      typedef typename make_const_reference<reference>::type const_reference;

      typedef typename reference_to_arg<reference>::type arg_type;
      typedef typename basic_type<reference>::type const_arg_type;

      explicit NoAliasProxy(reference Value) : Value_(Value) {}

      // compiler-synthesized copy-ctor is OK

      // is this also the copy-assignment operator?
      template <typename U>
      NoAliasProxy& operator=(U const& x)
      {
         assign(Value_, x);
         //Value_ = noalias(x);
         return *this;
      }

      // just in case...
      NoAliasProxy& operator=(NoAliasProxy const& x)
      {
         assign(Value_, x.value());
         //         Value_ = x;
         return *this;
      }

      // to allow assignment with temporary proxies,
      // computed assignment operators must be members here,
      // or see comment below for an alternative

      template <typename U>
      NoAliasProxy& operator+=(U const& x)
      {
         add(Value_, x);
         return *this;
      }

      template <typename U>
      NoAliasProxy& operator-=(U const& x)
      {
         subtract(Value_, x);
         return *this;
      }

      reference value() { return Value_; }
      const_reference value() const { return Value_; }

      operator NoAliasProxy<const_reference>() const
      { return NoAliasProxy<const_reference>(Value_); }


#if 0
      template <typename RangeType>
      typename VectorRange<const_arg_type, RangeType>::result_type
      range(RangeType r) const
      {
         using LinearAlgebra::range;
         return range(Value_, r);
      }

      template <typename RangeType>
      typename VectorRange<arg_type, RangeType>::result_type
      range(RangeType r)
      {
         using LinearAlgebra::range;
         return range(Value_, r);
      }

      template <typename IntType>
      typename boost::enable_if<exists<IntType>,
         VectorRange<const_arg_type, Range> >::type::result_type
      range(IntType first, IntType last) const
      {
         using LinearAlgebra::range;
         return range(Value_, Range(first, last));
      }

      template <typename IntType>
      typename boost::enable_if<exists<IntType>,
         VectorRange<arg_type, Range> >::type::result_type
      range(IntType first, IntType last)
      {
         using LinearAlgebra::range;
         return range(Value_, Range(first, last));
      }

      template <typename SliceType>
      typename VectorSlice<const_arg_type, SliceType>::result_type
      slice(SliceType s) const
      {
         using LinearAlgebra::slice;
         return slice(Value_, s);
      }

      template <typename SliceType>
      typename VectorSlice<arg_type, SliceType>::result_type
      slice(SliceType s)
      {
         using LinearAlgebra::slice;
         return slice(Value_, s);
      }

      template <typename StartType, typename SizeType, typename StrideType>
      typename boost::enable_if<exists<StartType>,
         VectorSlice<const_arg_type, Slice> >::type::result_type
      slice(StartType Start, SizeType Size, StrideType Stride) const
      {
         using LinearAlgebra::slice;
         return slice(Value_, Slice(Start, Size, Stride));
      }

      template <typename StartType, typename SizeType, typename StrideType>
      typename boost::enable_if<exists<StartType>,
         VectorSlice<arg_type, Slice> >::type::result_type
      slice(StartType Start, SizeType Size, StrideType Stride)
      {
         using LinearAlgebra::slice;
         return slice(Value_, Slice(Start, Size, Stride));
      }
#endif

   private:
      reference Value_;
};

// assignment etc

template <typename LHS, typename RHS>
struct Assign<NoAliasProxy<LHS&>&, RHS> : Assign<LHS&, RHS> {};

template <typename LHS, typename RHS>
struct Assign<LHS&, NoAliasProxy<RHS> >
   : Assign<LHS&, typename basic_type<RHS>::type> {};

template <typename LHS, typename RHS>
struct Assign<NoAliasProxy<LHS&>&, NoAliasProxy<RHS> >
   : Assign<LHS&, typename basic_type<RHS>::type> {};


template <typename LHS, typename RHS>
struct AssignCopy<NoAliasProxy<LHS&>&, RHS> : Assign<LHS&, RHS> {};

template <typename LHS, typename RHS>
struct AssignCopy<LHS&, NoAliasProxy<RHS> >
   : Assign<LHS&, typename basic_type<RHS>::type> {};

template <typename LHS, typename RHS>
struct AssignCopy<NoAliasProxy<LHS&>&, NoAliasProxy<RHS> >
   : Assign<LHS&, typename basic_type<RHS>::type> {};



template <typename LHS, typename RHS>
struct Add<NoAliasProxy<LHS&>&, RHS> : Add<LHS&, RHS> {};

template <typename LHS, typename RHS>
struct Add<LHS&, NoAliasProxy<RHS> >
   : Add<LHS&, typename basic_type<RHS>::type> {};

template <typename LHS, typename RHS>
struct Add<NoAliasProxy<LHS&>&, NoAliasProxy<RHS> >
   : Add<LHS&, typename basic_type<RHS>::type> {};


template <typename LHS, typename RHS>
struct AddCopy<NoAliasProxy<LHS&>&, RHS> : Add<LHS&, RHS> {};

template <typename LHS, typename RHS>
struct AddCopy<LHS&, NoAliasProxy<RHS> >
   : Add<LHS&, typename basic_type<RHS>::type> {};

template <typename LHS, typename RHS>
struct AddCopy<NoAliasProxy<LHS&>&, NoAliasProxy<RHS> >
   : Add<LHS&, typename basic_type<RHS>::type> {};



template <typename LHS, typename RHS>
struct Subtract<NoAliasProxy<LHS&>&, RHS> : Subtract<LHS&, RHS> {};

template <typename LHS, typename RHS>
struct Subtract<LHS&, NoAliasProxy<RHS> >
   : Subtract<LHS&, typename basic_type<RHS>::type> {};

template <typename LHS, typename RHS>
struct Subtract<NoAliasProxy<LHS&>&, NoAliasProxy<RHS> >
   : Subtract<LHS&, typename basic_type<RHS>::type> {};


template <typename LHS, typename RHS>
struct SubtractCopy<NoAliasProxy<LHS&>&, RHS> : Subtract<LHS&, RHS> {};

template <typename LHS, typename RHS>
struct SubtractCopy<LHS&, NoAliasProxy<RHS> >
   : Subtract<LHS&, typename basic_type<RHS>::type> {};

template <typename LHS, typename RHS>
struct SubtractCopy<NoAliasProxy<LHS&>&, NoAliasProxy<RHS> >
   : Subtract<LHS&, typename basic_type<RHS>::type> {};


// strip_noalias

template <typename T>
inline
T&
strip_noalias(NoAliasProxy<T>& x)
{
   return x.value();
}

template <typename T>
inline
T const&
strip_noalias(NoAliasProxy<T> const& x)
{
   return x.value();
}

//
// noalias function
//

template <typename T>
inline
typename boost::disable_if<is_proxy_reference<T>, NoAliasProxy<T&> >::type
noalias(T& x)
{
   return NoAliasProxy<T&>(x);
}

template <typename T>
inline
typename boost::disable_if<is_proxy_reference<T>, NoAliasProxy<T const&> >::type
noalias(T const& x)
{
   return NoAliasProxy<T const&>(x);
}

template <typename T>
inline
typename boost::enable_if<is_proxy<T>, NoAliasProxy<T> >::type
noalias(T x)
{
   return NoAliasProxy<T>(x);
}

template <typename T>
inline
NoAliasProxy<T>& noalias(NoAliasProxy<T>& x)
{
   return x;
}

template <typename T>
inline
NoAliasProxy<T> noalias(NoAliasProxy<T> x)
{
   return x;
}

#if 0
// this is the alternative scheme for handling temporary proxy references:
// we allow pass by value for the case where T is a proxy type.
template <typename T, typename U>
inline
NoAliasProxy<T>& operator+=(NoAliasProxy<T>& x, U const& y)
{
   x.value() += noalias(y);
   return x;
}

template <typename T, typename U>
inline
typename boost::enable_if<is_proxy<T>, NoAliasProxy<T> >::type
operator+=(NoAliasProxy<T> x, U const& y)
{
   x.value() += noalias(y);
}
#endif

} // namespace LinearAlgebra

#endif
