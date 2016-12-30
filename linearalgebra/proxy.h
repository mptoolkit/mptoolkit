// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/proxy.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

  proxy.h

  Some metafunctions to handle proxy types.
*/

#if !defined(PROXY_H_JDHCJKSDHTUI34756837478O)
#define PROXY_H_JDHCJKSDHTUI34756837478O

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/not.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/mpl/print.hpp>

namespace LinearAlgebra
{

//
// is_mutable_proxy
// metafunction returns true if T is a mutable proxy type
//
// To declare a type as a mutable proxy, either specialize is_mutable_proxy
// OR define a nested "typedef mpl::true_ T::proxy".
//
// a constant proxy declared as 'proxy<T const>' is mostly equivalent
// to 'proxy<T> const'.  That is, constness follows the logical const'ness
// of the referent.  But it is not possible in C++98 to distinguish a temporary
// 'proxy<T>' from a 'proxy<T> const', hence some functions that would
// accept a mutable proxy choose to perform a const-cast in this case.
// Thus, a mutable proxy can choose to specialize make_const_reference<T>
// to return 'proxy<T const>', rather than the default 'proxy<T> const'.
//
// C++11: now that we have r-value references, proxy<T const> is not necessary,
// we can just use proxy<T> const everywhere.
//

template <typename T>
struct is_mutable_proxy;

//
// is_const_proxy
// metafunction returns true if T is a constant proxy type
//
// To declare a type as a constant proxy, either specialize is_const_proxy
// OR define a nested "typedef mpl::true_ T::const_proxy".
//

template <typename T>
struct is_const_proxy;

//
// is_proxy
// metafunction returns true if T is a proxy type.
// This is simply is_mutable_proxy<T> || is_const_proxy<T>
//

template <typename T>
struct is_proxy;

//
// is_proxy_reference
// metafunction returns true if T is a proxy-reference,
// ie. either a proxy or a reference type.
// It should never be necessary to specialize this function.
//
// A reference to a proxy is usually erroneous and not a valid proxy-reference.
//

template <typename T>
struct is_proxy_reference;

template <typename T>
struct is_const_proxy_reference;

template <typename T>
struct is_mutable_proxy_reference;

//
// is_immediate
// metafunction returns true if T is not a proxy-reference,
// but is a 'small' value-type.  These objects are cheap to copy
// so we can treat them in some cases as similar to const proxies.
//
// TODO: possibly, this should be the default and we want to mark some
// types as proxyable?
//

template <typename T>
struct is_immediate;

template <typename T>
struct is_proxy_reference_or_immediate
   : boost::mpl::or_<is_proxy_reference<T>, is_immediate<T> > {};

template <typename T>
struct is_const_proxy_reference_or_immediate
   : boost::mpl::or_<is_const_proxy_reference<T>, is_immediate<T> > {};

//
// make_reference
// returns the (proxy) reference type.
// When constructing proxy objects, you want to
// store value types by reference, but proxy types
// are stored by the proxy itself (NOT a reference to the proxy -
// it is probably a temporary).
// It should never be necessary to specialize this function.
//

template <typename T>
struct make_reference;

//
// make_const_reference
//

template <typename T>
struct make_const_reference;

//
// proxy_lvalue
//
// if T is a mutable proxy, return T&,  Otherwise return T.
// This is a hack to work around not being able to bind
// a temporary proxy to a non-const reference.
//
// The parameter counterpart is proxy_lvalue_cast(T const& t),
// which does a const-cast if T is a proxy, and
// identity if T is not a proxy.
//

template <typename T, typename Enable = void>
struct proxy_lvalue
{
   typedef T type;
};

template <typename T>
struct proxy_lvalue<T, typename boost::enable_if<is_mutable_proxy<T> >::type>
{
   typedef T& type;
};

template <typename T, typename Enable = void>
struct proxy_lvalue_reference
{
   typedef T const& type;
};

template <typename T>
struct proxy_lvalue_reference<T, typename boost::enable_if<is_mutable_proxy<T> >::type>
{
   typedef T& type;
};

template <typename T>
inline
typename proxy_lvalue_reference<T>::type
proxy_lvalue_cast(T const& x)
{
   return const_cast<T&>(x);
}

// utility function - is something like this already in mpl?
template <typename T>
struct exists : boost::mpl::true_ { };

//
// implementation follows
//

namespace Private
{

template <typename T, typename Enable = void>
struct is_declared_proxy : boost::mpl::false_ {};

template <typename T>
struct is_declared_proxy<T, typename boost::enable_if<typename T::proxy>::type>
   : boost::mpl::true_
{
};

template <typename T, typename Enable = void>
struct is_declared_const_proxy : boost::mpl::false_ {};

template <typename T>
struct is_declared_const_proxy<T, typename boost::enable_if<typename T::const_proxy>::type>
   : boost::mpl::true_
{
};

template <typename T, typename Enable = void>
struct is_declared_immediate : boost::mpl::false_ {};

template <typename T>
struct is_declared_immediate<T, typename boost::enable_if<typename T::immediate>::type>
   : boost::mpl::true_
{
};

} // namespace Private

template <typename T>
struct is_mutable_proxy : Private::is_declared_proxy<T> {};

template <typename T>
struct is_mutable_proxy<T const> : boost::mpl::false_ {};

template <typename T>
struct is_const_proxy : Private::is_declared_const_proxy<T> {};

template <typename T>
struct is_const_proxy<T const>
   : boost::mpl::or_<Private::is_declared_proxy<T>,
                     Private::is_declared_const_proxy<T> > {};

template <typename T>
struct is_proxy : boost::mpl::or_<is_mutable_proxy<T>, is_const_proxy<T> > {};


template <typename T>
struct is_mutable_proxy_reference : is_mutable_proxy<T> {};

template <typename T>
struct is_mutable_proxy_reference<T&>
   : boost::mpl::not_<boost::mpl::or_<is_proxy<T>, boost::is_const<T> > > {};


template <typename T>
struct is_const_proxy_reference : is_const_proxy<T> {};

template <typename T>
struct is_const_proxy_reference<T&>
   : boost::mpl::and_<boost::mpl::not_<is_proxy<T> >, boost::is_const<T> > {};


template <typename T>
struct is_proxy_reference : is_proxy<T> {};

template <typename T>
struct is_proxy_reference<T&> : boost::mpl::not_<is_proxy<T> > {};


namespace Private
{
// make_const_ref_helper does the work of constructing the const proxy-reference
// type of T, which must not be a reference type nor have cv qualifiers.
template <typename T, typename Enable = void>
struct make_const_ref_helper  // primary declaration handles mutable proxies with no const_type
{
   typedef T const type;
};

// specialization for references, disable if is_proxy
template <typename T>
struct make_const_ref_helper<T, typename boost::disable_if<is_proxy<T> >::type>
{
   typedef T const& type;
};

// specialization for const proxies
template <typename T>
struct make_const_ref_helper<T, typename boost::enable_if<is_const_proxy<T> >::type>
{
   typedef T type;
};

// specialization for mutable proxies with a nested const_type
template <typename T>
struct make_const_ref_helper<T, typename boost::enable_if<
   boost::mpl::and_<
      is_mutable_proxy<T>
    , exists<typename T::const_type>
   >
>::type>
{
   typedef typename T::const_type type;
};

} // namespace Private


template <typename T>
struct make_const_reference : Private::make_const_ref_helper<T> {};

template <typename T>
struct make_const_reference<T&> : make_const_reference<T> {};

template <typename T>
struct make_const_reference<T const> : make_const_reference<T> {};


template <typename T>
struct make_reference : boost::mpl::if_<is_proxy<T>, T, T&> {};

template <typename T>
struct make_reference<T const> : make_const_reference<T> {};

template <typename T>
struct make_reference<T&> : make_reference<T> {};

} // namespace LinearAlgebra

#endif
