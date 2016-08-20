// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/hash_set.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

namespace ext
{

template <typename Key, class HashFun, class Cmp>
hash_set<Key, HashFun, Cmp>::hash_set()
  : TBase()
{
}

template <typename Key, class HashFun, class Cmp>
hash_set<Key, HashFun, Cmp>::hash_set(size_type n, hasher hf_, key_equal KEq)
  : TBase(n, hf_, KEq)
{
}

template <typename Key, class HashFun, class Cmp>
template <class Iter>
hash_set<Key, HashFun, Cmp>::hash_set(Iter first, Iter last)
  : TBase(first, last)
{
}

template <typename Key, class HashFun, class Cmp>
template <class Iter>
hash_set<Key, HashFun, Cmp>::hash_set(Iter first, Iter last, size_type n, hasher hf_, key_equal KEq)
  : TBase(first, last, n, hf_, KEq)
{
}

#if defined(USE_PSTREAM)

template <int Format, typename Key, class HashFun, class Cmp>
PStream::opstreambuf<Format>& operator<<(PStream::opstreambuf<Format>& out,
                                         hash_set<Key, HashFun, Cmp> const& hset)
{
   typename PStream::opstreambuf<Format>::size_type len = hset.size();
   out << len;
   for (typename hash_set<Key, HashFun, Cmp>::const_iterator I = hset.begin(); I != hset.end(); ++I)
   {
      out << *I;
   }
   return out;
}

template <int Format, typename Key, class HashFun, class Cmp>
PStream::ipstreambuf<Format>& operator>>(PStream::ipstreambuf<Format>& in,
                                         hash_set<Key, HashFun, Cmp>& hset)
{
   typename hash_set<Key, HashFun, Cmp>::value_type v;
   hset.clear();
   typename PStream::opstreambuf<Format>::size_type len;
   in >> len;
   for (int i = 0; i < len; ++i)
   {
      in >> v;
      hset.insert(v);
   }
   return in;
}

#endif // if defined(USE_PSTREAM)

} // namespace ext
