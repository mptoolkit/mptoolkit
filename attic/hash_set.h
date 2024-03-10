// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/hash_set.h
//
// Copyright (C) 2004-2021 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

/*
  hash_map.h

  A hash_set implementation

  Created 2004-02-24 Ian McCulloch
*/

#if !defined(HASH_SET_ZSDJFHUIWYR83479Y89UFY8934U)
#define HASH_SET_ZSDJFHUIWYR83479Y89UFY8934U

#include "private/hashtable.h"
#if defined(USE_PSTREAM)
#include "pstream/pstream.h"
#endif

namespace ext
{

// for a set, the value_type and key_type are the same.
// The ExtractKey and ConstructValue functions are no-ops.
template <typename T>
struct identity_func : public std::unary_function<T, T>
{
   T& operator()(T& x) const { return x; }
   T const& operator()(T const& x) const { return x; }
};

template <typename Key, class HashFun = hash<Key>, class Cmp = std::equal_to<Key> >
class hash_set : public Private::hash_table<Key,
                                            Key,
                                            HashFun,
                                            Cmp,
                                            identity_func<Key>,
                                            identity_func<Key> >
{
   private:
      typedef Private::hash_table<Key, Key, HashFun, Cmp,
                                  identity_func<Key>,
                                  identity_func<Key> >       TBase;
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
      typedef typename TBase::key_equal       value_equal; // not defined for hash_table
      typedef typename TBase::iterator        iterator;
      typedef typename TBase::const_iterator  const_iterator;
      typedef typename TBase::iterator        local_iterator;
      typedef typename TBase::const_iterator  const_local_iterator;

      hash_set();

      hash_set(size_type n, hasher hf_ = hasher(), key_equal KEq = key_equal());

      template <class Iter>
      hash_set(Iter first, Iter last);

      template <class Iter>
      hash_set(Iter first, Iter last, size_type n, hasher hf_ = hasher(), key_equal KEq = key_equal());
};

#if defined(USE_PSTREAM)
template <int Format, typename Key, class HashFun, class Cmp>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>&,
                                         hash_set<Key, HashFun, Cmp> const&);

template <int Format, typename Key, class HashFun, class Cmp>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>&,
                                         hash_set<Key, HashFun, Cmp>&);
#endif

} // namespace ext

#include "hash_set.cc"

#endif
