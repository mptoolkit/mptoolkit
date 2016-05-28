// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/vectorhashiterator.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

/*
  vectorhashiterator.h

  An adaptor to construct a linear algebra hashed iterator from
  a standard-style hashed container.

  Created 2005-01-06 Ian McCulloch
*/

#if !defined(VECTORHASHITERATOR_H_JKCHWUITY438T97UY89PU908UWP)
#define VECTORHASHITERATOR_H_JKCHWUITY438T97UY89PU908UWP

#include "common/trace.h"
#include "vectorinterface.h"
#include "scalar.h"
#include <iterator>
#include <boost/mpl/if.hpp>

namespace LinearAlgebra
{

// Container should be a hash_map-like container with a value_type of
// std::pair<size_type, value_type>.  ValueType should be
// either value_type or value_type const.  If it is value_type const, 
// this behaves as a const-iterator, otherwise it is a non-const iterator.
template <typename Container, typename ValueType = typename Container::value_type::second_type>
class VectorHashIterator
{
   public:
      typedef Container container_type;
      typedef typename boost::mpl::if_<
	 boost::is_const<ValueType>, 
	 typename container_type::const_iterator,
	 typename container_type::iterator>::type container_iterator;
      typedef typename boost::mpl::if_<
	 boost::is_const<ValueType>, 
	 Container const,
	 Container>::type ContainerType;
      typedef ContainerType& container_reference;
      typedef ContainerType* container_pointer;

      typedef typename boost::remove_const<ValueType>::type value_type;
      typedef ValueType& reference;
      typedef ValueType* pointer;
      typedef vector_iterator_hashed category;

      VectorHashIterator() {}

      explicit VectorHashIterator(container_reference Data)
	 : Data_(&Data), Here_(Data.begin()) {}

      VectorHashIterator(container_reference Data, container_iterator const& Here)
	 : Data_(&Data), Here_(Here) {}

      // this ctor handles conversion from non-const to const.
      template <typename OtherVal>
      VectorHashIterator(VectorHashIterator<Container, OtherVal> const& x)
	 : Data_(&x.data()), Here_(x.base_iterator()) {}

      VectorHashIterator& operator++()
      { DEBUG_CHECK(Here_ != this->end()); ++Here_; return *this; }

      VectorHashIterator operator++(int)
      { DEBUG_CHECK(Here_ != this->end()); return VectorHashIterator(*Data_, Here_++); }

      VectorHashIterator& operator--()
      { DEBUG_CHECK(Here_ != Data_->begin()); --Here_; return *this; }

      VectorHashIterator operator--(int)
      { DEBUG_CHECK(Here_ != Data_->begin()); return VectorHashIterator(*Data_, Here_--); }

      size_type index() const { DEBUG_CHECK(Here_ != this->end()); return Here_->first; }

      reference operator*() const { return Here_->second; }

      pointer operator->() const { return &Here_->second; }

      operator bool() const { return Here_ != this->end(); }

      // members for the hashed_iterator concept
      size_type size() const { return Data_->size(); }
      reference operator()(size_type n) const;
      bool element_exists(size_type n) const;

      container_reference data() const { return *Data_; }
      container_iterator base_iterator() const { return Here_; }

   private:
      container_iterator end() const { return Data_->end(); }

      container_pointer Data_;
      container_iterator Here_;
};

template <typename Container, typename ValueType>
struct VectorHashIteratorHelper
{
   typedef typename VectorHashIterator<Container, ValueType>::container_pointer container_pointer;
   typedef typename VectorHashIterator<Container, ValueType>::reference reference;
   
   static reference lookup(container_pointer p, size_type n)
   {
      return (*p)[n];
   }
};

template <typename Container, typename ValueType>
struct VectorHashIteratorHelper<Container, ValueType const>
{
   typedef typename VectorHashIterator<Container, ValueType const>::value_type value_type;
   typedef typename VectorHashIterator<Container, ValueType const>::container_pointer container_pointer;
   typedef typename VectorHashIterator<Container, ValueType const>::container_iterator container_iterator;
   typedef typename VectorHashIterator<Container, ValueType const>::reference reference;
   
   static reference lookup(container_pointer p, size_type n)
   {
      container_iterator I = p->find(n);
      if (I == p->end())
      {
	 return static_zero_or_die<value_type>();
      }
      else return I->second;
   }
};

template <typename Container, typename ValueType>
inline
typename VectorHashIterator<Container, ValueType>::reference
VectorHashIterator<Container, ValueType>::operator()(size_type n) const
{
   return VectorHashIteratorHelper<Container, ValueType>::lookup(Data_, n);
}

template <typename Container, typename ValueType>
inline
bool
VectorHashIterator<Container, ValueType>::element_exists(size_type n) const
{
   return Data_->count(n);
}

} // namespace LinearAlgebra

#endif
