// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/hash_map.h
//
// Copyright (C) 2002-2021 Ian McCulloch <ian@qusim.net>
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
  hash_map.h

  A simple hash_map implementation

  Created 2002-11-15 Ian McCulloch

*/

#if !defined(HASH_MAP_SDJFHUIWYR83479Y89UFY8934U)
#define HASH_MAP_SDJFHUIWYR83479Y89UFY8934U

#include "private/hashtable.h"
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

namespace ext
{

// select1st is potentially useful in its own right
template <typename T>
struct select1st : public std::unary_function<T, typename T::first_type>
{
   typename T::first_type& operator()(T& x) const { return x.first; }
   typename T::first_type const& operator()(T const& x) const { return x.first; }
};

// make1st constructs a pair from a copy of the first_type and
// a default-constructed second-type.
template <typename T>
struct make1st : public std::unary_function<typename T::first_type, T>
{
   T operator()(typename T::first_type const& x) { return T(x, typename T::second_type()); }
};

template <typename Key,
          typename T,
          class HashFun = hash<Key>,
          class Cmp = std::equal_to<Key> >
class hash_map : public Private::hash_table<Key,
                                            std::pair<Key, T>,
                                            HashFun,
                                            Cmp,
                                            select1st<std::pair<Key, T> >,
                                            make1st<std::pair<Key, T> > >
{
   private:
      typedef Private::hash_table<Key, std::pair<Key, T>, HashFun, Cmp,
                                  select1st<std::pair<Key, T> >,
                                  make1st<std::pair<Key, T> >  >       TBase;
   public:
      typedef typename TBase::key_type        key_type;
      typedef typename TBase::value_type      value_type;
      typedef typename TBase::pointer         pointer;
      typedef typename TBase::const_pointer   const_pointer;
      typedef typename TBase::reference       reference;
      typedef typename TBase::const_reference const_reference;
      typedef typename TBase::size_type       size_type;
      typedef typename TBase::difference_type difference_type;
      typedef typename TBase::hasher          hasher;
      typedef typename TBase::key_equal       key_equal;
      typedef typename TBase::iterator        iterator;
      typedef typename TBase::const_iterator  const_iterator;
      typedef typename TBase::iterator        local_iterator;
      typedef typename TBase::const_iterator  const_local_iterator;
      typedef T                               mapped_type;

      hash_map();

      hash_map(size_type n, hasher hf_ = hasher(), key_equal KEq = key_equal());

      template <class Iter>
      hash_map(Iter first, Iter last);

      template <class Iter>
      hash_map(Iter first, Iter last, size_type n,
               hasher hf_ = hasher(), key_equal KEq = key_equal());

      mapped_type& operator[](key_type const& k);
};

#if defined(USE_PSTREAM)
template <int Format, typename Key, typename T, class HashFun, class Cmp>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>&,
                                         hash_map<Key, T, HashFun, Cmp> const&);

template <int Format, typename Key, typename T, class HashFun, class Cmp>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>&,
                                         hash_map<Key, T, HashFun, Cmp>&);
#endif

} // namespace ext

#include "hash_map.cc"

#endif
