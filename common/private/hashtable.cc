// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/private/hashtable.cc
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

#include <math.h>

namespace ext
{

namespace Private
{

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::hash_table()
  : hf(hasher()), KeyCompare(key_equal()), KeyExtractor(EKey()), ValueCtor(CValue()),
    Size(0), RehashSize(size_t(Private::DefaultMaxLoadFactor * Private::MinBuckets)),
    NumBuckets(Private::MinBuckets), MaxLoadFactor(Private::DefaultMaxLoadFactor)
#if defined(HASHTABLE_POWER2)
    , Mask(Private::MinMask)
#endif
{
   BucketVector.assign(NumBuckets, BucketType(NodeList.end(), 0));
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::hash_table(size_type n,
                                                           hasher hf_,
                                                           key_equal KEq,
                                                           EKey KExtract,
                                                           CValue VCtor)
  : hf(hf_), KeyCompare(KEq), KeyExtractor(KExtract), ValueCtor(VCtor), Size(0),
    NumBuckets(Private::FindNextPrime(n)),
    MaxLoadFactor(Private::DefaultMaxLoadFactor)
#if defined(HASHTABLE_POWER2)
    , Mask(Private::MinMask)
#endif
{
   RehashSize = Private::DefaultMaxLoadFactor * NumBuckets;
   BucketVector.assign(n, BucketType(NodeList.end(), 0));
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
template <class Iter>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::hash_table(Iter first, Iter last)
  : hf(hasher()), KeyCompare(key_equal()), KeyExtractor(EKey()), ValueCtor(CValue()),
    Size(0), RehashSize(size_type(Private::DefaultMaxLoadFactor * Private::MinBuckets)),
  NumBuckets(Private::MinBuckets), MaxLoadFactor(Private::DefaultMaxLoadFactor)
#if defined(HASHTABLE_POWER2)
    , Mask(Private::MinMask)
#endif
{
   BucketVector.assign(NumBuckets, BucketType(NodeList.end(), 0));
   while (first != last)
   {
      insert(*first++);
   }
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
template <class Iter>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::hash_table(Iter first,
                                                           Iter last,
                                                           size_type n,
                                                           hasher hf_,
                                                           key_equal KEq,
                                                           EKey KExtract,
                                                           CValue VCtor)
  : hf(hf_), KeyCompare(KEq), KeyExtractor(KExtract), ValueCtor(VCtor), Size(0),
    NumBuckets(Private::FindNextPrime(n)),
    MaxLoadFactor(Private::DefaultMaxLoadFactor)
#if defined(HASHTABLE_POWER2)
    , Mask(Private::MinMask)
#endif
{
   RehashSize = Private::DefaultMaxLoadFactor * NumBuckets;
   BucketVector.assign(n, BucketType(NodeList.end(), 0));
   while (first != last)
   {
      insert(*first++);
   }
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::clear()
{
   NodeList.clear();
   std::fill(BucketVector.begin(), BucketVector.end(), BucketType(NodeList.end(), 0));
   Size = 0;
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
std::pair<typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::iterator, bool>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::insert(value_type const& v)
{
   BucketType* Bucket = &BucketVector[bucket(KeyExtractor(v))];
   iterator Where = this->FindInBucket(Bucket, KeyExtractor(v));
   if (Where == NodeList.end())
   {
      if (Size >= RehashSize)
      {
         this->Grow();
         Bucket = &BucketVector[bucket(KeyExtractor(v))];
      }

      Where = Bucket->first = NodeList.insert(Bucket->first, v);
      ++Bucket->second;
      ++Size;
      return std::make_pair(Where, true);
   }
   else
   {
      return std::make_pair(Where, false);
   }
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
template <class Iter>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::insert(Iter first, Iter last)
{
   while (first != last)
   {
      this->insert(*first);
      ++first;
   }
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::erase(key_type const& k)
{
   BucketType* Bucket = &BucketVector[this->bucket(k)];
   iterator Where = this->FindInBucket(Bucket, k);

   if (Where == NodeList.end()) return;  // not found
   DoErase(Bucket, Where);
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
inline
typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::iterator
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::find(key_type const& k)
{
   return this->FindInBucket(&BucketVector[this->bucket(k)], k);
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
inline
typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::const_iterator
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::find(key_type const& k) const
{
   return this->FindInBucket(&BucketVector[this->bucket(k)], k);
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::reference
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::find_or_insert(key_type const& k)
{
   BucketType* Bucket = &BucketVector[bucket(k)];
   iterator Where = this->FindInBucket(Bucket, k);

   if (Where == NodeList.end())
   {
      if (Size >= RehashSize)
      {
         this->Grow();
         Bucket = &BucketVector[bucket(k)];
      }

      Where = Bucket->first = NodeList.insert(Bucket->first, ValueCtor(k));
      ++Bucket->second;
      ++Size;
   }
   return *Where;
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
std::pair<typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::iterator,
          typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::iterator>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::equal_range(key_type const& k)
{
   iterator Where = this->find(k);
   iterator WhereEnd = Where;
   if (WhereEnd != NodeList.end()) ++WhereEnd;
   return std::make_pair(Where, WhereEnd);
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
std::pair<typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::const_iterator,
          typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::const_iterator>
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::equal_range(key_type const& k) const
{
   const_iterator Where = this->find(k);
   const_iterator WhereEnd = Where;
   if (WhereEnd != NodeList.end()) ++WhereEnd;
   return std::make_pair(Where, WhereEnd);
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::rehash(size_type NewBuckets)
{
#if defined(HASHTABLE_POWER2)
   // round NewBuckets up to a power of 2, and determine the Mask
   Mask = 0;
   size_t TBuckets = NewBuckets-1;
   NewBuckets = 1;
   while (TBuckets != 0)
   {
      TBuckets >>= 1;
      NewBuckets <<= 1;
      Mask <<= 1;
      Mask |= 1;
   }
#else
   NewBuckets = Private::FindNextPrime(NewBuckets);
#endif

   // Early return if the new bucket count is the same as the old
   if (NewBuckets == NumBuckets) return;

   // if we are going to throw at all, it had better be in the next line...
   BucketVectorType(NewBuckets, BucketType(NodeList.end(), size_type(0))).swap(BucketVector);

   // from now on we don't throw

   NumBuckets = NewBuckets;
   RehashSize = size_t(NumBuckets * MaxLoadFactor);

   // Rearrange the NodeList and fill the BucketVector.
   // The strategy to do this is to take the front element in NodeList
   // and move it to the correct position, and repeat Size times.
   DEBUG_CHECK(NodeList.size() == Size);
   for (int i = 0; i < int(Size); ++i)
   {
      iterator I = NodeList.begin();
      BucketType* Bucket = &BucketVector[bucket(KeyExtractor(*I))];
      NodeList.splice(Bucket->first, NodeList, I);
      Bucket->first = I;
      ++Bucket->second;
   }
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::reserve(size_type n)
{
   size_type DesiredBuckets = size_type(ceil(n / MaxLoadFactor));
   if (DesiredBuckets > NumBuckets) rehash(DesiredBuckets);
}

//
// helper functions
//

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::iterator
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::FindInBucket(BucketType* Bucket, key_type const& k)
{
   size_type Count = Bucket->second;
   if (Count == 0) // early return if the node is empty
   {
      return NodeList.end();
   }
   // else
   iterator I = Bucket->first;
   if (KeyCompare(k, KeyExtractor(*I))) return I;
   ++I;
   --Count;
   while (Count != 0 && !(KeyCompare(k, KeyExtractor(*I))))
   {
      ++I;
      --Count;
   }
   return (Count != 0) ? I : NodeList.end();
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
typename hash_table<Key, T, HashFun, Cmp, EKey, CValue>::const_iterator
hash_table<Key, T, HashFun, Cmp, EKey, CValue>::FindInBucket(BucketType const* Bucket, key_type const& k) const
{
   size_type Count = Bucket->second;
   if (Count == 0) // early return if the node is empty
   {
      return NodeList.end();
   }
   // else
   const_iterator I = Bucket->first;
   if (KeyCompare(k, KeyExtractor(*I))) return I;
   ++I;
   --Count;
   while (Count != 0 && !(KeyCompare(k, KeyExtractor(*I))))
   {
      ++I;
      --Count;
   }
   return (Count != 0) ? I : NodeList.end();
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::DoErase(BucketType* Bucket, iterator Where)
{
   // There is a choice as to what iterator we use for empty buckets.
   // If we use NodeList.end() then we don't need to do anything special when inserting
   // an item into a previously empty bucket, but we have to reset the iterator to
   // NodeList.end() if we remove the last item from a bucket.
   // The alternative is to have empty bucket iterators point at nowhere in particular,
   // which makes erase() slightly faster at the expense of making insert() slightly slower.
   // We choose the former case.

   --Bucket->second;
   --Size;

   if (Bucket->second == 0)
   {
      NodeList.erase(Where);
      Bucket->first = NodeList.end();
   }
   else if (Where == Bucket->first)
   {
      Bucket->first = NodeList.erase(Where);
   }
   else
   {
      NodeList.erase(Where);
   }
}

template <typename Key, typename T, class HashFun, class Cmp, class EKey, class CValue>
void hash_table<Key, T, HashFun, Cmp, EKey, CValue>::Grow()
{
   this->rehash(size_t(NumBuckets * 1.5));  // the choice of primes means the actual growth rate is ~ 2
}

} // namespace Private

} // namespace ext
