// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/index.h
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
  index.h

  type Index for indirect views of containers.

  Created 2004-04-13 Ian McCulloch
*/

#if !defined(INDEX_H_CHSJKH4UIWHUIHNLUIREHT8793Y8794YWFE80)
#define INDEX_H_CHSJKH4UIWHUIHNLUIREHT8793Y8794YWFE80

#include "vectoroperationsbase.h"
#include "noalias.h"
#include "crtp_vector.h"

namespace LinearAlgebra
{

//
// IndirectIterator
//
// Iterator over a vector that is indexed by another iterator.
// The base iterator must be a dense iterator.  The indexing iterator
// can be anything.
//
// We don't have any debug checks here; the underlying iterators should do that
// for us.
//

template <typename BaseIter, typename IndexIterType>
class IndirectIterator
{
   public:
      typedef BaseIter                                base_iterator;
      typedef IndexIterType                           index_iterator;

      BOOST_MPL_ASSERT((boost::is_integral<typename index_iterator::value_type>));
      BOOST_MPL_ASSERT((boost::is_convertible<typename base_iterator::category, vector_iterator_dense>));

      typedef typename base_iterator::value_type          value_type;
      typedef typename base_iterator::reference           reference;
      typedef typename base_iterator::pointer             pointer;
      typedef typename index_iterator::category category;

      IndirectIterator() {}

      template <typename B>
      IndirectIterator(B const& Base, index_iterator const& IndexIter)
	: Base_(Base), IndexIter_(IndexIter) {}

      template <typename Other, typename OtherIndex>
      IndirectIterator(IndirectIterator<Other, OtherIndex> const& Other_)
	 : Base_(Other_.base()), IndexIter_(Other_.index_iter()) {}

      template <typename Other, typename OtherIndex>
      IndirectIterator& operator=(IndirectIterator<Other, OtherIndex> const& Other_)
      { Base_ = Other_.base(); IndexIter_ = index_iter(); return *this; }

      IndirectIterator& operator++() { ++IndexIter_; return *this; }
      IndirectIterator operator++(int) { return IndirectIterator(Base_, IndexIter_++); }

      IndirectIterator& operator+=(difference_type n) { IndexIter_ += n; return *this; }

      reference operator*() const { return Base_[*IndexIter_]; }

      reference operator[](difference_type n) const { return Base_[IndexIter_[n]]; }

      pointer operator->() const 
	 //         { base_iterator b(Base_); b += *IndexIter_; return b.operator->(); }
         { return &Base_[*IndexIter_]; }

      reference operator()(difference_type n) const { return Base_[IndexIter_(n)]; }

      size_type index() const { return IndexIter_.index(); }

      operator bool() const { return bool(IndexIter_); }

      base_iterator& base() { return Base_; }
      base_iterator const& base() const { return Base_; }

      index_iterator const& index_iter() const { return IndexIter_; }

   //      size_type loc() const { return Loc_; }

   private:
      base_iterator Base_;
      index_iterator IndexIter_;
};

// we never want to instantiate a IndirectIterator<IndirectIterator<T> >,
// better to compose the functors together instead.
// so make it an error to do so
// 2006-03-03: well, maybe it isn't really that evil after all.  Without it,
// we have to implement specializations for all functions that produce
// a IndirectVector.
#if 0
template <typename BaseIterType, typename F1, typename F2>
struct IndirectIterator<IndirectIterator<BaseIterType, F2> , F1>
{
   //   BOOST_MPL_ASSERT((boost::mpl::false_));
};
#endif

template <typename VecType, typename IndexType>
class IndirectVector;

template <typename VecType, typename IndexType>
class IndirectVector : public VectorBase<IndirectVector<VecType, IndexType> >
{
   public:
      BOOST_MPL_ASSERT((is_proxy_reference<VecType>));
      BOOST_MPL_ASSERT((is_const_proxy_reference<IndexType>));

      typedef is_mutable_proxy_reference<VecType> proxy;
      typedef is_const_proxy_reference<VecType> const_proxy;

      typedef VecType vector_reference;
      typedef IndexType index_type;

      typedef typename make_const_reference<index_type>::type  const_index_reference;
      typedef typename make_const_reference<vector_reference>::type  const_vector_reference;

      IndirectVector(vector_reference Vec, const_index_reference IndexVec)
	 : Vec_(Vec), IndexVec_(IndexVec) 
      {
         using LinearAlgebra::size;
	 DEBUG_PRECONDITION(size(IndexVec) == 0U || min(IndexVec) >= 0); 
	 DEBUG_PRECONDITION(size(IndexVec) == 0U || size_type(max(IndexVec)) < size(Vec));
      }

      template <typename OtherVec, typename OtherIndex>
      IndirectVector(IndirectVector<OtherVec, OtherIndex> const& Other)
	 : Vec_(Other.base()), IndexVec_(Other.index_vector())
      {
      }

      IndirectVector& operator=(IndirectVector const& x)
      {
         using LinearAlgebra::size;
	 typedef typename interface<VecType>::value_type value_type;
	 //std::vector<value_type> Temp(size(x));
         // noalias(Temp) = x;
	 assign(*this, x);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, IndirectVector&>::type 
      operator=(U const& x)
      {
         using LinearAlgebra::size;
	 typedef typename interface<VecType>::value_type value_type;
	 // TODO: the vector type should be something better here!
	 std::vector<value_type> Temp(size(x));
	 noalias(Temp) = x;

	 assign(*this, Temp);
	 return *this;
      }

      template <typename U>
      typename boost::enable_if<is_vector<U>, IndirectVector&>::type 
      operator=(NoAliasProxy<U> const& x)
      {
	 assign(*this, x.value());
	 return *this;
      }
   
      vector_reference base() { return Vec_; }
      const_vector_reference base() const { return Vec_; }
      const_index_reference index_vector() const { return IndexVec_; }

   private:
      vector_reference Vec_;
      index_type IndexVec_;
      size_type Size_;
};

// interface

template <typename T, typename I, typename Tv>
struct IndirectVectorInterface {};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, VECTOR_EXPRESSION(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef VECTOR_EXPRESSION(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, LOCAL_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef LOCAL_VECTOR(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, COMPRESSED_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef COMPRESSED_VECTOR(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, INJECTIVE_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef INJECTIVE_VECTOR(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, HASHED_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef HASHED_VECTOR(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, ORDERED_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef ORDERED_VECTOR(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, SEARCH_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef SEARCH_VECTOR(Tv, T) type;
};

template <typename T, typename Vv, typename Vi, typename Tv>
struct IndirectVectorInterface<T, DENSE_VECTOR(Vv, Vi), Tv>
{
   typedef Tv value_type;
   typedef DENSE_VECTOR(Tv, T) type;
};

// Cutoff at Dense - we cannot be a stride or contiguous vector
// with an arbitary indirect index.

template <typename VecType, typename IndexType>
struct interface<IndirectVector<VecType, IndexType> > 
   : IndirectVectorInterface<IndirectVector<VecType, IndexType>, 
			     typename interface<IndexType>::type,
			     typename interface<VecType>::value_type> {};

// iterators

template <typename VecType, typename IndexType>
struct Iterate<IndirectVector<VecType, IndexType>&>
{
   typedef IndirectIterator<typename iterator<VecType>::type, 
			    typename const_iterator<IndexType>::type> result_type;
   typedef IndirectVector<VecType, IndexType>& argument_type;
   result_type operator()(argument_type x) const
   {
      return result_type(iterate(x.base()), iterate(x.index_vector()));
   }
};

template <typename VecType, typename IndexType>
struct Iterate<IndirectVector<VecType, IndexType> >
{
   typedef IndirectIterator<typename const_iterator<VecType>::type, 
			    typename const_iterator<IndexType>::type> result_type;
   typedef IndirectVector<VecType, IndexType> const& argument_type;
   typedef typename basic_type<IndexType>::type BasicIndexType;
   result_type operator()(argument_type x) const
   {
      return result_type(Iterate<VecType>()(x.base()), Iterate<BasicIndexType>()(x.index_vector()));
   }
};

// size

template <typename VecType, typename IndexType>
struct Size<IndirectVector<VecType, IndexType> >
{
   typedef size_type result_type;
   typedef IndirectVector<VecType, IndexType> const& argument_type;
   result_type operator()(argument_type v) const { return size(v.index_vector()); }
};

// get_element

template <typename VecType, typename IndexType>
struct GetVectorElement<
   IndirectVector<VecType, IndexType>
 , typename boost::enable_if<
      is_defined<
         GetVectorElement<typename basic_type<IndexType>::type>
      > 
   >::type
>
{
   typedef typename GetVectorElement<typename basic_type<VecType>::type>::result_type result_type;
   typedef IndirectVector<VecType, IndexType> const& first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(first_argument_type v, second_argument_type i) const
   { return get_element(v.base(), get_element(v.index_vector(), i)); }
};

template <typename VecType, typename IndexType>
struct GetVectorElement<
   IndirectVector<VecType, IndexType>&
 , typename boost::enable_if<
      boost::mpl::and_<
         is_defined<GetVectorElement<typename basic_type<IndexType>::type> >
       , is_defined<GetVectorElement<typename basic_type<VecType>::type&> >
      > 
   >::type
>
{
   typedef typename GetVectorElement<typename basic_type<VecType>::type&>::result_type result_type;
   typedef IndirectVector<VecType, IndexType>& first_argument_type;
   typedef size_type second_argument_type;
   result_type operator()(first_argument_type v, second_argument_type i) const
   { return get_element(v.base(), get_element(v.index_vector(), i)); }
};

// TODO: set_element, add_element, subtract_element, multiply_element

template <typename T, typename U, 
	  typename TInterface = typename interface<T>::type,
	  typename UInterface = typename interface<U>::type>
struct Index {};

template <typename T, typename U>
inline
typename boost::enable_if<typename boost::mpl::and_<is_mutable_proxy_reference<T>,
						    is_defined<Index<T&, U> > >,
			  Index<T&, U> >::type::result_type
index(T const& x, U const& y)
{
   return Index<T&, U>()(const_cast<T&>(x),y);
}

template <typename T, typename U>
inline
typename boost::disable_if<typename boost::mpl::and_<is_mutable_proxy_reference<T>,
						     is_defined<Index<T&, U> > >,
			  Index<T, U> >::type::result_type
index(T const& x, U const& y)
{
   return Index<T, U>()(x,y);
}

template <typename T, typename U>
inline
typename Index<T&, U>::result_type
index(T& x, U const& y)
{
   return Index<T&, U>()(x,y);
}

// default implementation.
template <typename T, typename U, typename IntType, typename Enable = void>
struct ImplementIndex {};

template <typename T, typename U, typename IntType>
struct ImplementIndex<T, U, IntType, typename boost::enable_if<boost::is_integral<IntType> >::type>
{
   typedef IndirectVector<T, U> result_type;
   typedef T first_argument_type;
   typedef U second_argument_type;
   result_type operator()(T x, U y) const { return result_type(x,y); }
};

template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct Index<T, U, DENSE_VECTOR(Tv, Ti), LOCAL_VECTOR(Uv, Ui)>
   : ImplementIndex<typename make_const_reference<T>::type, 
		    typename make_const_reference<U>::type, Uv> {};

template <typename T, typename U, typename Tv, typename Ti, typename Uv, typename Ui>
struct Index<T&, U, DENSE_VECTOR(Tv, Ti), LOCAL_VECTOR(Uv, Ui)>
   : ImplementIndex<typename make_reference<T>::type, 
		    typename make_const_reference<U>::type, Uv> {};



// VectorBracket overload

template <typename T, typename Arg, 
          typename Tv, typename Ti, 
          typename Uv, typename Ui>
struct VectorBracketInterface<T, Arg, DENSE_VECTOR(Tv, Ti), LOCAL_VECTOR(Uv, Ui)>
   : Index<T, Arg> {};



} // namespace LinearAlgebra

#endif
