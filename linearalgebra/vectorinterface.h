// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectorinterface.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  vectorinterface.h

  New version of the vector interface heirachy

  Created 2004-12-24 Ian McCulloch
*/

#if !defined(VECTORINTERFACE_H_DHJFKHURWEIRY4572893475489YRUI34897)
#define VECTORINTERFACE_H_DHJFKHURWEIRY4572893475489YRUI34897

#include "common/trace.h"
#include "interface.h"
#include "vectorfwd.h"
#include "operations.h"
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cstddef>
#include <vector>

namespace LinearAlgebra
{

//
// The vector interface heirachy
//
// VectorExpression
// | |
// | Distributed
// |
// Local                          - defines iterators
// | |
// | Compressed                   - generic sparse vector
// | |
// | Injective                    - sparse vector, unique indices
// | | |
// | | Hashed                     - hashed vector, O(1) lookup
// | |
// | Ordered                      - sparse vector, ordered by index
// | |
// | Search                       - ordered vector with O(log N) lookup
// |
// Dense                          - generic dense vector, random access
// |
// Stride                         - direct pointer access, but non-trivial stride
// |
// Contiguous                     - direct pointer access, contiguous
//

// Arguably, Dense should inherit from Search, but
// that would imply sparse algorithms would always be useable with
// dense matrices/vectors. But that isn't necessarily desirable; for example
// assigning sparse = dense shouldn't be allowed (need a filter function
// to eliminate zeros), but it would be hard to avoid with that scheme.

// VectorExpression defines only size(), eval() and assignment semantics.
// Local defines the Iterate function, and iterator types.
//
// Currently, it is not possible for a vector to satisfy the Local      
// interface without also satisfying either the Dense or Compressed      
// interface.  
// 
// The interface satisfied by all iterators is an indexed iterator,
// with dereference, a member function index() that returns the 
// index of the iterator location, and a conversion to bool
// that returns true if the iterator is non-singular, false otherwise.
// Whether this is a forward iterator or bidirectional is yet to be determined;
// with the current implementation, not all derived iterators allow bidirectional.
//
// In addition, dense vectors define a dense iterator, that
// specifies the indices start from 0 and are contiguous,
// and allows random access operator[] with an *offset* relative to the
// current iterator location.
//
// Hashed sparse iterators define average-case O(1) random access
// via operator(), which uses absolute indexing and is not relative
// to the current iterator location.  (the use of () instead of []
// for the direct access reflects this change in semantics.)
//
// Usage of the macros is recommended over directly naming the interface
// types; the interface types are verbose nested typedefs and using
// them directly prohibits later modifications to the interface heirachy.
// The exception is possibly VectorExpression, which is always at top level.
//
// Note that the iterator heirachy (defined below) is somewhat different
// to the interface heirachy.  In particular, dense iterators are a
// refinement of (sparse) ordered iterators.
//
// The distributed vector interface is templated on an interface
// that represents the components stored locally.
//
// TODO: how do fixed-size vectors fit into this heirachy? An extra parameter to the
// dense interface?
//
// TODO: we probably need another type for sparse vectors with a fixed structure.
//
// TODO: an is_injective function would be useful.
//

namespace Concepts
{

// we write the names in full here so we can use templates with commas in them
// as arguments to the TEMPLATEX macros inside the interface macros.

#define TEMPLATE2(T,U) T,U /**/
#define TEMPLATE3(T,U,V) T,U,V /**/
#define TEMPLATE4(T,U,V,W) T,U,V,W /**/
#define TEMPLATE5(T,U,V,W,X) T,U,V,W,X /**/

template <typename Value, typename T>
struct VectorExpression
{
   typedef Value value_type;
   typedef T     interface_type;
};

#define ANY_VECTOR(Value, T) LinearAlgebra::Concepts::VectorExpression<Value, T > /**/

#define VECTOR_EXPRESSION(Value, T) LinearAlgebra::Concepts::VectorExpression<Value, T > /**/

template <typename LocalPart, typename T>
struct DistributedVector {};

#define DISTRIBUTED_VECTOR(Value, Local, T)				\
  LinearAlgebra::Concepts::VectorExpression< Value, 			\
   LinearAlgebra::Concepts::DistributedVector< Local, T > > /**/


template <typename T>
struct Local {};

#define LOCAL_STORAGE(T) 			\
   LinearAlgebra::Concepts::Local< T > /**/

#define LOCAL_VECTOR(Value, T)				\
LinearAlgebra::Concepts::VectorExpression< Value ,	\
   LinearAlgebra::Concepts::Local< T > > /**/


template <typename T>
struct Compressed {}; // : Local<Compressed<T> > {};

#define COMPRESSED_STORAGE(T) 				\
   LinearAlgebra::Concepts::Local<			\
      LinearAlgebra::Concepts::Compressed< T > > /**/

#define COMPRESSED_VECTOR(Value, T)			\
LinearAlgebra::Concepts::VectorExpression< Value ,	\
   LinearAlgebra::Concepts::Local<			\
      LinearAlgebra::Concepts::Compressed< T > > > /**/


template <typename T>
struct Injective {}; // : Compressed<Injective<T> > {};

#define INJECTIVE_STORAGE(T) 					\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::Compressed<			\
         LinearAlgebra::Concepts::Injective< T > > > /**/

#define INJECTIVE_VECTOR(Value, T)				\
LinearAlgebra::Concepts::VectorExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::Compressed<			\
         LinearAlgebra::Concepts::Injective< T > > > > /**/


template <typename T>
struct Hashed {}; // : Injective<Hashed<T> > {};

#define HASHED_STORAGE(T)					\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::Compressed<			\
         LinearAlgebra::Concepts::Injective<			\
            LinearAlgebra::Concepts::Hashed< T > > > > /**/

#define HASHED_VECTOR(Value, T)					\
LinearAlgebra::Concepts::VectorExpression< Value ,		\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::Compressed<			\
         LinearAlgebra::Concepts::Injective<			\
            LinearAlgebra::Concepts::Hashed< T > > > > > /**/


template <typename T>
struct Ordered {}; // : Injective<Ordered<T> > {};

#define ORDERED_STORAGE(T)					\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::Compressed<			\
         LinearAlgebra::Concepts::Injective<			\
            LinearAlgebra::Concepts::Ordered< T > > > > /**/

#define ORDERED_VECTOR(Value, T)					\
LinearAlgebra::Concepts::VectorExpression< Value ,			\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::Compressed<				\
         LinearAlgebra::Concepts::Injective<				\
            LinearAlgebra::Concepts::Ordered< T > > > > > /**/


template <typename T>
struct Search {}; // : Ordered<Search<T> > {};

#define SEARCH_STORAGE(T)						\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::Compressed<				\
         LinearAlgebra::Concepts::Injective<				\
            LinearAlgebra::Concepts::Ordered<				\
               LinearAlgebra::Concepts::Search< T > > > > > /**/

#define SEARCH_VECTOR(Value, T)						\
LinearAlgebra::Concepts::VectorExpression< Value ,			\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::Compressed<				\
         LinearAlgebra::Concepts::Injective<				\
            LinearAlgebra::Concepts::Ordered<				\
               LinearAlgebra::Concepts::Search< T > > > > > > /**/


template <typename T>
struct Dense {}; // : Local<Dense<T> > {};

#define DENSE_STORAGE(T) 				\
   LinearAlgebra::Concepts::Local<			\
      LinearAlgebra::Concepts::Dense< T > > /**/

#define DENSE_VECTOR(Value, T)				\
LinearAlgebra::Concepts::VectorExpression< Value ,	\
   LinearAlgebra::Concepts::Local<			\
      LinearAlgebra::Concepts::Dense< T > > > /**/


template <typename T>
struct Stride {}; // : Dense<Stride<T> > {};

#define STRIDE_STORAGE(T) 				\
   LinearAlgebra::Concepts::Local<			\
      LinearAlgebra::Concepts::Dense<			\
         LinearAlgebra::Concepts::Stride< T > > > /**/

#define STRIDE_VECTOR(Value, T)				\
LinearAlgebra::Concepts::VectorExpression< Value ,	\
   LinearAlgebra::Concepts::Local<			\
      LinearAlgebra::Concepts::Dense<			\
         LinearAlgebra::Concepts::Stride< T > > > > /**/


template <typename T>
struct Contiguous {}; // : Stride<Contiguous<T> > {};

#define CONTIGUOUS_STORAGE(T)					\
   LinearAlgebra::Concepts::Local<				\
      LinearAlgebra::Concepts::Dense<				\
         LinearAlgebra::Concepts::Stride<			\
            LinearAlgebra::Concepts::Contiguous< T > > > > /**/

#define CONTIGUOUS_VECTOR(Value, T)					\
LinearAlgebra::Concepts::VectorExpression< Value ,			\
   LinearAlgebra::Concepts::Local<					\
      LinearAlgebra::Concepts::Dense<					\
         LinearAlgebra::Concepts::Stride<				\
            LinearAlgebra::Concepts::Contiguous< T > > > > > /**/

} // namespace Concepts

// RebindInterface

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, VECTOR_EXPRESSION(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef VECTOR_EXPRESSION(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, LOCAL_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef LOCAL_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, COMPRESSED_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef COMPRESSED_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, INJECTIVE_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef INJECTIVE_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, HASHED_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef HASHED_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, ORDERED_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef ORDERED_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, SEARCH_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef SEARCH_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, DENSE_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef DENSE_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, STRIDE_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef STRIDE_VECTOR(Tv, void) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct RebindInterface<T, CONTIGUOUS_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef CONTIGUOUS_VECTOR(Tv, void) type;
};

//
// is_vector, boolean function to determine if a type
// has a Local interface.
//

namespace Private
{

template <typename T>
struct is_vector_helper : boost::mpl::false_ {};

template <typename S, typename T>
struct is_vector_helper<VECTOR_EXPRESSION(S,T)> : boost::mpl::true_ {};

} // namespace Private

template <typename T, typename Enable = void>
struct is_vector : boost::mpl::false_
{
};

template <typename T>
struct is_vector<T, typename boost::enable_if<exists<
   typename interface<typename boost::remove_const<T>::type>::type> >::type>
   : Private::is_vector_helper<
   typename interface<typename boost::remove_const<T>::type>::type> {};

//
// when allocating temporaries, it is useful to know
// simply whether a type should be sparse or dense.
//

struct vector_abstract_dense 
{
   typedef vector_abstract_dense type;
};

struct vector_abstract_sparse 
{
   typedef vector_abstract_sparse type;
};

// union of vectors; dense+anything -> dense;  sparse+sparse -> sparse.
template <typename T, typename U>
struct vector_abstract_or 
   : boost::mpl::if_<boost::is_same<T, U>, T, vector_abstract_dense> {};

// intersection vectors; sparse+anything -> sparse;  dense+dense -> dense.
template <typename T, typename U>
struct vector_abstract_and
   : boost::mpl::if_<boost::is_same<T, U>, T, vector_abstract_sparse> {};

//
// make_vector_from_abstract
//
// provides a 'default' choice of a concrete vector type based on
// the abstract type (ie. either sparse or dense).
//

template <typename T, typename Abstract>
struct make_vector_from_abstract;

template <typename T>
struct make_vector_from_abstract<T, vector_abstract_dense>
{
   typedef Vector<T> type;
};

template <typename T>
struct make_vector_from_abstract<T, vector_abstract_sparse>
{
   typedef MapVector<T> type;
};

template <typename T, typename Sv, typename Si>
struct abstract_interface_interface<T, DENSE_VECTOR(Sv, Si)>
   : vector_abstract_dense {};

template <typename T, typename Sv, typename Si>
struct abstract_interface_interface<T, COMPRESSED_VECTOR(Sv, Si)>
   : vector_abstract_sparse {};

// make_value_from_interface for vectors

template <typename T, typename Sv, typename Si>
struct make_value_from_interface<T, ANY_VECTOR(Sv, Si)>
   : make_vector_from_abstract<Sv, typename abstract_interface<T>::type> {};

//
// iterators
//

// FIXME: iterator and const_iterator should work as 'expected' for all combinations
// of top-level const and references.

template <typename T>
struct Iterate;

template <typename T, typename Enable = void>
struct iterator
{
   typedef typename Iterate<T>::result_type type;
   typedef typename Iterate<T>::result_type const_type;
};

#if 0
template <typename T>
struct iterator<T const>
{
   typedef typename Iterate<T>::result_type type;
   typedef typename Iterate<T>::result_type const_type;
};
#endif

template <typename T>
struct iterator<T, typename boost::enable_if<
   exists<typename Iterate<T&>::result_type> >::type>
{
   typedef typename Iterate<T&>::result_type type;
   typedef typename Iterate<T>::result_type const_type;
};

template <typename T>
struct iterator<T&, typename boost::enable_if<
   exists<typename Iterate<T&>::result_type> >::type>
{
   typedef typename Iterate<T&>::result_type type;
   typedef typename Iterate<T>::result_type const_type;
};

template <typename T, typename Enable = void>
struct const_iterator
{
   typedef typename Iterate<typename basic_type<T>::type>::result_type type;
};

#if defined(USE_FORCE_INLINE)
template <typename T>
typename Iterate<T&>::result_type
iterate(T& x) __attribute__((always_inline));

template <typename T>
typename Iterate<typename proxy_lvalue<T>::type>::result_type
iterate(T const& x) __attribute__((always_inline));
#endif

template <typename T>
inline
typename Iterate<T&>::result_type
iterate(T& x)
{
   return Iterate<T&>()(x);
}

template <typename T>
inline
typename Iterate<typename proxy_lvalue<T>::type>::result_type
iterate(T const& x)
{
   return Iterate<typename proxy_lvalue<T>::type>()(proxy_lvalue_cast(x));
}

template <typename T>
inline
typename boost::disable_if<is_defined<Iterate<typename proxy_lvalue<T>::type> >,
                           Iterate<T> >::type::result_type
iterate(T const& x)
{
   return Iterate<T>()(x);
}

template <typename T>
inline
typename Iterate<T>::result_type
const_iterate(T const& x)
{
   return Iterate<T>()(x);
}

// iterator_traits are probably unnecessary; the
// iterators must be classes anyway.

template <typename T, typename Enable = void>
struct iterator_traits
{
   typedef typename T::value_type value_type;
   typedef typename T::reference  reference;
   typedef typename T::pointer    pointer;
   typedef typename T::category   category;
};

//
// iterator categories
//

// generic sparse unordered iterator
struct vector_iterator_sparse {};

// unordered but distinct indices
struct vector_iterator_injective : vector_iterator_sparse {};

// hashed indices; contains operator() for lookup, 
// member function has_element(), member function size().
// Size can only be an approximation in some cases; is it actually useful for anything?
struct vector_iterator_hashed : vector_iterator_injective {};

// sparse but ordered indices.  No size() member.
struct vector_iterator_ordered : vector_iterator_injective {};

// dense ordered indices (we don't allow for dense but unordered)
// Has operator[], member function size()
struct vector_iterator_dense : vector_iterator_ordered {};

} // namespace LinearAlgebra

#endif
