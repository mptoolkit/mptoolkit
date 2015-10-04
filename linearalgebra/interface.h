/* -*- C++ -*- $Id$

  interface.h

  top level declarations for defining interfaces and proxies.

  Created 2004-12-20 Ian McCulloch
*/

#if !defined(INTERFACE_H_JDHCJKSDHTUI34756837478O)
#define INTERFACE_H_JDHCJKSDHTUI34756837478O

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
#include "common/trace.h"
#include "proxy.h"

namespace LinearAlgebra
{

typedef std::size_t size_type;
typedef std::ptrdiff_t difference_type;

//
// interface - a function for getting the interface (concept) corresponding to a type.
// If the type is a container (vector or matrix), then it also contains a typedef
// value_type, which is the (concrete) value_type that it contains.
//

// AnyScalar: catch-all for the interface of a type that is not
// known to the LinearAlgebra library, so we regard it as a scalar
template <typename T> struct AnyScalar;

template <typename T, typename Enable = void>
struct interface
{
   //   typedef void type;
   typedef AnyScalar<T> type;
   typedef T value_type;
};

template <typename T>
struct interface<T const> : interface<T> {};

template <typename T>
struct interface<T&> : interface<T> {};

//
// make_value - a function for getting the value type of a const-reference
// or proxy reference.  Conceptually, it returns the basic type of
// an expression.  Generally, it shouldn't need to be specialized itself,
// it defaults to stripping of (const) references, and then 
// for proxies, returning make_value_from_interface<T>, or
// for non-proxies, T itself.
// 
// make_value_from_interface - generally T will be a proxy type here.
// This function generates a value_type from an interface type.
//

template <typename T, typename TInterface = typename interface<T>::type>
struct make_value_from_interface {};

template <typename T, typename U>
struct make_value_from_interface<T, AnyScalar<U>>
{
   typedef T type;  // or U?
};

template <typename T>
struct make_value 
   : boost::mpl::if_<is_proxy<T>, make_value_from_interface<T>, boost::mpl::identity<T> >::type {};

template <typename T>
struct make_value<T&> : make_value<T> {};

template <typename T>
struct make_value<T const> : make_value<T> {};

//
// make_value_with_zero
// make a value_type that also has a zero_value.  If type T already has a zero
// representation (specified by has_zero<T> - which is equivalent to 
// Zero<T>::result_type being defined, then make_value_with_zero<T> is equivalent
// to make_value<T>.  Otherwise, it has a nested type which is
// value_with_zero<make_value<T>::type>.
//

template <typename T, typename Enable = void>
struct make_value_with_zero;

// has_zero

template <typename T, typename Enable = void>
struct has_zero;

template <typename T, typename Enable = void>
struct Zero;

template <typename T, typename Enable = void>
struct IsZero;

//
// result_value
//
// metafunction to return the value_type corresponding to the result of a metafunction.
// A short-cut for make_value<T::result_type>.
// For example, 
// typename result_value<Negate<T> >::type x = negate(y);
// We cannot use Negate<T>::result_type directly here, as it may be a proxy.
// This is effectively a replacement for
// decltype(negate(x)), since we cannot use decltype with expression templates.
//

template <typename T>
struct result_value : make_value<typename T::result_type> {};

//
// call_type
//
// determines the type to use to pass parameters.
// inpired by boost/call_traits
//

template <typename T, typename Enable = void>
struct call_type
{
   typedef T const& type;
};

template <typename T>
struct call_type<T, typename boost::enable_if<is_proxy_reference<T> >::type>
{
   typedef T type;
};

//
// arg_type
//
// 'inverse' of call_type
//

template <typename T, typename Enable = void>
struct arg_type
{
   typedef T type;
};

template <typename T>
struct arg_type<T const&>
{
   typedef T type;
};

// reference_to_arg - convert a reference to an argument type suitable as a template
// argument.  Undefined if T is not a reference.
//
// T (mutable proxy) -> T& 
// T (const proxy) -> T
// T& -> T&
// T const& -> T

template <typename T, typename Enable = void>
struct reference_to_arg;

template <typename T>
struct reference_to_arg<T&>
{
   typedef T& type;
};

template <typename T>
struct reference_to_arg<T const&>
{
   typedef T type;
};

template <typename T>
struct reference_to_arg<T, typename boost::enable_if< is_mutable_proxy<T> >::type>
{
   typedef T& type;
};

template <typename T>
struct reference_to_arg<T, typename boost::enable_if< is_const_proxy<T> >::type>
{
   typedef T type;
};

template <typename T>
struct reference_to_arg<T const, typename boost::enable_if< is_const_proxy<T const> >::type>
{
   typedef T type;
};

// basic_type - is this the same thing as boost::decay ?

template <typename T, typename Enable = void>
struct basic_type
{
   typedef T type;
};

template <typename T>
struct basic_type<T&> : basic_type<T> {};

template <typename T>
struct basic_type<T const> : basic_type<T> {};


//
// value_type
//
// for T a vector or matrix, value_type<T>::type is the logical type of
// the container.  It is an alias for interface<T>::value_type.
// If T is not a vector nor matrix, then value_type<T> is an empty struct.
//

template <typename T, typename Enable = void>
struct value_type {};

template <typename T>
struct value_type<T, typename boost::enable_if<exists<typename interface<T>::value_type> >::type>
{
   typedef typename interface<T>::value_type type;
};


// RebindInterface - FIXME : is this used?

template <typename T, typename OldInterface, typename Tv = typename T::value_type>
struct RebindInterface 
{
   typedef void type;
};

//
// abstract_interface
// This is used by matrix and vector expressions to get a type of the final
// expression (eg for use as a temporary)
//

// abstract_interface_interface gets the abstract interface from the 'real'
// interface.
// The easiest way to define it is to declare nested typedef 
// T::abstract_type.  Otherwise, specialize vector_abstract_of_expression itself.

template <typename T, typename TInterface = typename interface<T>::type>
struct abstract_interface_interface
{
   typedef typename T::abstract_interface type;
};

template <typename T>
struct abstract_interface : abstract_interface_interface<T, typename interface<T>::type> {};

// For convenience, strip const references (since this is often called on proxy types)
template <typename T>
struct abstract_interface<T&> : abstract_interface<T> {};

template <typename T>
struct abstract_interface<T const> : abstract_interface<T> {};

//
// implies
// helper metafunction - doesn't really belong here

template <typename T, typename U>
struct implies : boost::mpl::or_<boost::mpl::not_<T>, boost::mpl::and_<T,U> > {};

} // namespace LinearAlgebra

#endif
