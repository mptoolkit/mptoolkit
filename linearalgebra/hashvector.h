// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/hashvector.h
//
// Copyright (C) 2005-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  hashvector.h

  A sparse vector based on ext::hash_map

  Created 2005-01-07 Ian McCulloch

  We could improve on the implementation by avoiding the
  requirement for a default constructor in some places.
*/

#if !defined(HASHVECTOR_H_UIYFUIH57YT793W8Y9)
#define HASHVECTOR_H_UIYFUIH57YT793W8Y9

#include "vectorhashiterator.h"
#include "vectoroperations.h"
#include "common/hash_map.h"
#include "vectoraddition.h"
#include "vectortransform.h"
#include "crtp_vector.h"

namespace LinearAlgebra
{

template <typename T>
class HashVector : public VectorBase<HashVector<T> >
{
   private:
      typedef ext::hash_map<size_type, T> container_type;
      typedef typename container_type::iterator base_iterator;
      typedef typename container_type::const_iterator base_const_iterator;

   public:
      typedef T value_type;
      typedef VectorHashIterator<container_type, value_type> iterator;
      typedef VectorHashIterator<container_type, value_type const> const_iterator;
      typedef T& reference;
      typedef T const& const_reference;

      HashVector() : Size_(0) {}

      explicit HashVector(size_type Size) : Size_(Size) {}

      template <typename U>
      HashVector(U const& x, typename boost::enable_if<is_vector<U> >::type* dummy = 0);

      template <typename U>
      HashVector(NoAliasProxy<U> const& x, typename boost::enable_if<is_vector<U> >::type* dummy = 0);

      template <typename U>
      typename boost::enable_if<is_vector<U>, HashVector<T>&>::type
      operator=(U const& x);

      template <typename U>
      typename boost::enable_if<is_vector<U>, HashVector<T>&>::type
      operator=(NoAliasProxy<U> const& x);

      const_reference operator[](size_type n) const;
      reference operator[](size_type n);

      size_type size() const { return Size_; }

      const_iterator iterate() const { return const_iterator(Data_); }
      const_iterator c_iterate() const { return const_iterator(Data_); }
      iterator iterate() { return iterator(Data_); }

      iterator iterate_at(size_type n) 
      { return iterator(Data_, Data_.find(n)); }

      const_iterator iterate_at(size_type n) const 
      { return const_iterator(Data_, Data_.find(n)); }

      size_type nnz() const { return Data_.size(); }

      template <typename U>
      void set_element(size_type n, U const& x);

      template <typename U>
      void add_element(size_type n, U const& x);

      template <typename U>
      void subtract_element(size_type n, U const& x);

      void zero_element(size_type n);

      void zero_all() { Data_.clear(); }

      void resize(size_type NewSize) { Data_.clear(); Size_ = NewSize; }

      void swap(HashVector& Other)
      { std::swap(Size_, Other.Size_); Data_.swap(Other.Data_); }

   private:
      size_type Size_;
      container_type Data_;
};

// interface

template <typename T>
struct interface<HashVector<T> >
{
   //   typedef HASHED_VECTOR(T, HashVector<T>) type;
   typedef HASHED_VECTOR(T, void) type;
   typedef T value_type;
};

// iterators

template <typename T>
struct Iterate<HashVector<T>&>
{
   typedef HashVector<T>& argument_type;
   typedef typename HashVector<T>::iterator result_type;
   result_type operator()(argument_type x) const
   {
      return x.iterate();
   }
};

template <typename T>
struct Iterate<HashVector<T> >
{
   typedef HashVector<T> const& argument_type;
   typedef typename HashVector<T>::const_iterator result_type;
   result_type operator()(argument_type x) const
   {
      return x.iterate();
   }
};

// resize

template <typename T>
struct Resize<HashVector<T>&>
{
   typedef void result_type;
   typedef HashVector<T>& first_argument_type;
   typedef size_type second_argument_type;
   void operator()(HashVector<T>& v, size_type n) const
   {
      v.resize(n);
   }
};

// nnz

template <typename T>
struct NNZ<MapVector<T> >
{
   typedef HashVector<T> const& argument_type;
   typedef size_type result_type;

   result_type operator()(argument_type x) const
   {
      return x.nnz();
   }
};

// operator+=, operator-=

template <typename Scalar, typename U>
inline
HashVector<Scalar>&
operator+=(HashVector<Scalar>& x, U const& y)
{
   // TODO: the choice of HashVector as temporary type isn't necessarily optimal here
   HashVector<Scalar> Temp(y);
   add(x, Temp);
   return x;
}

template <typename Scalar, typename U>
inline
HashVector<Scalar>&
operator+=(HashVector<Scalar>& x, NoAliasProxy<U> const& y)
{
   add(x, y.value());
   return x;
}

template <typename Scalar, typename U>
inline
HashVector<Scalar>&
operator-=(HashVector<Scalar>& x, U const& y)
{
   // TODO: the choice of HashVector as temporary type isn't necessarily optimal here
   HashVector<Scalar> Temp(y);
   subtract(x, Temp);
   return x;
}

template <typename Scalar, typename U>
inline
HashVector<Scalar>&
operator-=(HashVector<Scalar>& x, NoAliasProxy<U> const& y)
{
   subtract(x, y.value());
   return x;
}

} // namespace LinearAlgebra

#include "hashvector.cc"

#endif
