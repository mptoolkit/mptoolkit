// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/private/hashtable.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ian@qusim.net>
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
  hashtable.h

  Generic hash-table interface, as a common base of hash_map and hash_set.

  Converted from the initial hash_map implementation on 2004-02-23 Ian McCulloch

  hash_map version created 2002-11-15 Ian McCulloch

  This implementation uses a doubly-linked list as the underlying container,
  with a std::vector<std::pair<iterator, size_t> > to denote the buckets, with
  iterator pointing to the first element in the bucket, and size_t elements
  in the bucket.

  Overview of implementation:

  the BucketVector becomes a vector of begin/end pairs of iterators into the NodeList

  Initially all pairs point to NodeList.end()

  If we add an item to a bucket which is empty, then the bucket moves to the
  front of NodeList:
    last = NodeList.begin();
    NodeList.push_front(NewValue);
    first = NodeList.begin();
    BucketVector[bucket] = make_pair(first, last);

  If we add an item to a non-empty bucket, we must add it at the back of the bucket.

  The preceeding two rules guarantee that elements in the same bucket stay adjacent.

  Erasing elements works fine, except we need to ensure that if the bucket becomes empty
  we reset the begin/end pair to be NodeList.end().  Otherwise our iterator might become
  invalidated.

  PROBLEM: what happens if we erase the beginning of a bucket() ???

  Solution:
    Instead of keeping a begin/end pair for each bucket, we keep only the begin, and a count.
    that way there is no end iterator to be invalidated, and it is safe to insert before
    the previous begin.
    In fact, in this scheme it is safe to make new buckets at the end of the NodeList.

    But, the requirements say that end(n) must be constant-time.  How?

    Solution: Keep the buckets ordered.  Cannot do this without scanning?

    Alternate solution: keep a pair (first, last-1) for each bucket.  Empty buckets are
    treated specially.

  The proposed standard suggests that it is possible to call erase() with a const_iterator.
  However it seems to be not possible to erase() a const_iterator in a list, therefore
  we require a non-const iterator.

*/

#if !defined(HASHTABLE_H_HIUHIUF743Y8732T7EGGDG3287)
#define HASHTABLE_H_HIUHIUF743Y8732T7EGGDG3287

#include <functional>
#include <vector>
#include <list>
#include <algorithm>
#include "common/trace.h"
#include "common/hash.h"

namespace ext
{

namespace Private
{

// KeyT is the key type that is used in compare operations.
// ValueT is the actual type that is stored in the container.
// HashFun is a hash function for objects of type KeyT.
// KeyCmp is a comparison function for objects of type KeyT.
// ExtractKey is a functor that returns an object of type KeyT
// from an object of type ValueT.
// ConstructValue is a functor that constructs a ValueT given
// a KeyT.  This needs to set a policy for any components
// of ValueT that cannot be determined from the KeyT.  The only
// purpose of this parameter is to make operator[] of hash_map
// slightly more efficient (otherwise find_or_insert() would
// have to take a value_type as a parameter instead of a key_type,
// meaning that every call to operator[] would require constructing
// a ValueT even if it is subsequently not used).

template <typename KeyT, typename ValueT, class HashFun,
          class KeyCmp, class ExtractKey, class ConstructValue>
class hash_table
{
   public:
      typedef KeyT                    key_type;
      typedef ValueT                  value_type;
      typedef value_type*             pointer;
      typedef const value_type*       const_pointer;
      typedef value_type&             reference;
      typedef const value_type&       const_reference;
      typedef size_t                  size_type;
      typedef ptrdiff_t               difference_type;
      typedef HashFun                 hasher;
      typedef KeyCmp                  key_equal;

   private:
      typedef std::list<value_type>   NodeListType;

   public:
      typedef typename NodeListType::iterator iterator;
      typedef typename NodeListType::const_iterator const_iterator;

      typedef typename NodeListType::iterator local_iterator;
      typedef typename NodeListType::const_iterator const_local_iterator;

      hash_table();

      hash_table(size_type n, hasher hf_ = hasher(), key_equal KEq = key_equal(),
                 ExtractKey KExtractor = ExtractKey(), ConstructValue VCtor = ConstructValue());

      template <class Iter>
      hash_table(Iter first, Iter last);

      template <class Iter>
      hash_table(Iter first, Iter last, size_type n, hasher hf_ = hasher(),
                 key_equal KEq = key_equal(), ExtractKey KExtractor = ExtractKey(),
                 ConstructValue VCtor = ConstructValue());

      iterator begin() { return NodeList.begin(); };
      iterator end() { return NodeList.end(); }

      const_iterator begin() const { return NodeList.begin(); };
      const_iterator end() const { return NodeList.end(); }

      local_iterator begin(size_type n) { return BucketVector[n].first; }
      local_iterator end(size_type n)
      { local_iterator Temp = BucketVector[n].first; advance(Temp, BucketVector[n].second); return Temp; }

      const_local_iterator begin(size_type n) const { return BucketVector[n].first; }
      const_local_iterator end(size_type n) const
      { const_local_iterator Temp = BucketVector[n].first; advance(Temp, BucketVector[n].second); return Temp; }

      size_type size() const { return Size; }
      size_type max_size() const { return BucketVector.max_size(); }
      size_type count(key_type const& k) const { return this->find(k) == this->end() ? 0 : 1; }

      size_type bucket_count() const { return NumBuckets; }
      size_type max_bucket_count() const { return BucketVector.back().max_size(); }

#if defined(HASHTABLE_POWER2)
      size_type bucket(key_type const& k) const { return hf(k) & Mask; }
#else
      size_type bucket(key_type const& k) const { return hf(k) % NumBuckets; }
#endif

      size_type bucket_size(size_type b) const { return BucketVector[b].second; }

      double load_factor() const { return Size / BucketVector.size(); }
      double max_load_factor() const { return MaxLoadFactor; }
      void set_max_load_factor(double f) { MaxLoadFactor = f; RehashSize = size_t(MaxLoadFactor * NumBuckets); }

      hasher hash_funct() const { return hf; }
      key_equal key_eq() const { return KeyCompare; }

      bool empty() const { return Size == 0; }

      void clear();

      std::pair<iterator, bool> insert(value_type const& v);

      template <class Iter>
      void insert(Iter first, Iter last);

      void erase(key_type const& k);
      void erase(iterator I) { DoErase(&BucketVector[bucket(KeyExtractor(*I))], I); }

      iterator find(key_type const& k);
      const_iterator find(key_type const& k) const;

      reference find_or_insert(key_type const& k);

      std::pair<iterator, iterator> equal_range(key_type const& k);
      std::pair<const_iterator, const_iterator> equal_range(key_type const& k) const;

      // rehashes with a new bucket count of at least n
      void rehash(size_type n);

      // adds buckets such that there is room for at least n elements
      // in the container before another rehash.
      // Does not ever reduce the bucket count.
      void reserve(size_type n);

   //      size_type collisions() const { return Collisions; }

   private:
      typedef typename NodeListType::iterator NodeListIterType;

      typedef std::pair<NodeListIterType, size_type> BucketType;
      typedef NodeListIterType BucketIterType;

      typedef std::vector<BucketType> BucketVectorType;

      typedef typename BucketVectorType::iterator BucketVectorIterType;

      // returns an iterator to the location of k in Bucket, or end() if no such key
      iterator FindInBucket(BucketType* Bucket, key_type const& k);
      const_iterator FindInBucket(BucketType const* Bucket, key_type const& k) const;

      void DoErase(BucketType* Bucket, iterator Where);

      void Grow();

      hasher hf;
      key_equal KeyCompare;
      ExtractKey KeyExtractor;
      ConstructValue ValueCtor;

      size_type Size;       // number of items in the container
      size_type RehashSize; // once we hit this size, it is time to rehash
      size_type NumBuckets;       // number of buckets.

      double MaxLoadFactor;

      NodeListType NodeList;
      BucketVectorType BucketVector;

#if defined(HASHTABLE_POWER2)
      size_t Mask;
#endif
};

} // namespace Private

} // namespace ext

#include "hashtable.cc"

#endif
