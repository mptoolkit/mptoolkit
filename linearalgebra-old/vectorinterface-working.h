// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectorinterface-working.h
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
// | DistributedVector
// |
// LocalVector                    - defines iterators
// | |
// | CompressedVector             - generic sparse vector
// | |
// | InjectiveVector              - sparse vector, unique indices
// | | |
// | | HashedVector               - hashed vector, O(1) lookup
// | |
// | OrderedVector                - sparse vector, ordered by index
// | |
// | SearchVector                 - ordered vector with O(log N) lookup
// |
// DenseVector                    - generic dense vector, random access
// |
// StrideVector                   - direct pointer access, but non-zero stride
// |
// ContiguousVector               - direct pointer access, contiguous
//

// Arguably, DenseVector should inherit from SearchVector, but
// that would imply sparse matrix algorithms would always be useable with
// dense matrices. But that isn't necessarily desirable; for example
// assigning sparse = dense shouldn't be allowed (need a filter function
// to eliminate zeros), but it would be hard to avoid with that scheme.

// VectorExpression defines only size() and eval() semantics.
// TODO: the relationship between assign and eval is a bit confused.
// LocalVector defines the Iterate function, and iterator types.
//
// Currently, it is not possible for a vector to satisfy the LocalVector
// interface without also satisfying either the DenseVector or CompressedVector
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
// TODO: how do fixed-size vectors fit into this heirachy?
//
// TODO: The intention was that DistributedVector would be templated on
// an interface type that represents the part of the vector which is stored locally.
//
// TODO: we probably need another type for sparse vectors with a fixed structure.
//

// this indirect stuff is an experiment to get
// the derived-most type as a typedef of the interface type without
// causing an instantiation loop.

template <typename T>
struct indirect
{
   typedef T type;
};

template <typename Value, typename T>
struct VectorExpression
{
   typedef Value value_type;
   typedef T     interface_type;
};

template <typename Value, typename T>
struct indirect<VectorExpression<Value, T> >
{
   typedef typename indirect<T>::type type;
};

#define ANY_VECTOR(Value, T) VectorExpression<Value, T > /**/
#define VECTOR_EXPRESSION(Value, T) VectorExpression<Value, T > /**/

template <typename Value, typename T>
struct DistributedVector : public VectorExpression<Value, DistributedVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

template <typename Value, typename T>
struct indirect<DistributedVector<Value, T> >
{
   typedef typename indirect<T>::type type;
};

#define DISTRIBUTED_VECTOR(Value, T) VectorExpression< Value, DistributedVector< Value, T > > /**/

template <typename Value, typename T>
struct LocalVector : public VectorExpression<Value, LocalVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define LOCAL_VECTOR(Value, T) VectorExpression< Value, LocalVector< Value, T > > /**/

template <typename Value, typename T>
struct CompressedVector : public LocalVector<Value, CompressedVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define COMPRESSED_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, CompressedVector< Value, T > > > /**/

template <typename Value, typename T>
struct InjectiveVector : public CompressedVector<Value, InjectiveVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define INJECTIVE_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, CompressedVector< Value, \
 InjectiveVector< Value, T > > > > /**/

template <typename Value, typename T>
struct HashedVector : public InjectiveVector< Value, HashedVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define HASHED_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, CompressedVector< Value, \
 InjectiveVector< Value, HashedVector< Value, T > > > > > /**/

template <typename Value, typename T>
struct OrderedVector : public InjectiveVector<Value, OrderedVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define ORDERED_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, CompressedVector< Value, \
 InjectiveVector< Value, OrderedVector< Value, T > > > > > /**/

template <typename Value, typename T>
struct SearchVector : public OrderedVector<Value, SearchVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define SEARCH_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, CompressedVector< Value, \
 InjectiveVector< Value, OrderedVector< Value, SearchVector< Value, T > > > > > > /**/

template <typename Value, typename T>
struct DenseVector : public LocalVector<Value, DenseVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define DENSE_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, DenseVector< Value, T > > > /**/

template <typename Value, typename T>
struct StrideVector : public DenseVector<Value, StrideVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define STRIDE_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, DenseVector< Value, \
 StrideVector< Value, T > > > > /**/

template <typename Value, typename T>
struct ContiguousVector : public StrideVector<Value, ContiguousVector<Value, T> >
{
   //   typedef typename indirect<T>::type type;
   //   typedef typename T::type type;
};

#define CONTIGUOUS_VECTOR(Value, T) \
 VectorExpression< Value, LocalVector< Value, DenseVector< Value, \
 StrideVector< Value, ContiguousVector< Value, T > > > > > /**/

//
// is_vector, boolean function to determine if a type
// has a LocalVector interface.
//

namespace Private
{

template <typename T>
struct is_vector_helper : public boost::mpl::false_ {};

template <typename S, typename T>
struct is_vector_helper<VectorExpression<S, T> > : public boost::mpl::true_ {};

} // namespace Private

template <typename T, typename Enable = void>
struct is_vector : public boost::mpl::false_
{
};

template <typename T>
struct is_vector<T, typename boost::enable_if<exists<
   typename interface<typename boost::remove_const<T>::type>::type> >::type>
   : public Private::is_vector_helper<
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

//
// vector_abstract_of_expression
//
// For a type that is an expression, return the abstract type (sparse or dense).
// The easiest way to define it is to declare nested typedef
// T::abstract_type.  Otherwise, specialize vector_abstract_of_expression itself.
//

template <typename T, typename TInterface = typename interface<T>::type>
struct abstract_interface
{
   typedef typename T::abstract_interface type;
};

template <typename T, typename Sv, typename Si>
struct abstract_interface<T, DENSE_VECTOR(Sv, Si)>
   : vector_abstract_dense {};

template <typename T, typename Sv, typename Si>
struct abstract_interface<T, COMPRESSED_VECTOR(Sv, Si)>
   : vector_abstract_dense {};

// make_value_from_interface for vectors

template <typename T, typename Sv, typename Si>
struct make_value_from_interface<T, ANY_VECTOR(Sv, Si)>
   : make_vector_from_abstract<Sv, typename abstract_interface<T>::type> {};

//
// iterators
//

template <typename T>
struct Iterate;

template <typename T, typename Enable = void>
struct iterator
{
   typedef typename Iterate<T const>::result_type type;
   typedef typename Iterate<T const>::result_type const_type;
};

template <typename T>
struct iterator<T, typename boost::enable_if<
   exists<typename Iterate<T>::result_type> >::type>
{
   typedef typename Iterate<T>::result_type type;
   typedef typename Iterate<T const>::result_type const_type;
};

template <typename T, typename Enable = void>
struct const_iterator
{
   typedef typename Iterate<T const>::result_type type;
};

template <typename T>
inline
typename Iterate<T>::result_type
iterate(T& x)
{
   return Iterate<T>()(x);
}

template <typename T>
inline
typename Iterate<T const>::result_type
iterate(T const& x)
{
   return Iterate<T const>()(x);
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
