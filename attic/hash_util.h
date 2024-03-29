// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/hash_util.h
//
// Copyright (C) 2002-2021 Ian McCulloch <ian@qusim.net>
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
  hash_util.h

  Some functions to compliment hash_map

  Created 2002-11-17 Ian McCulloch
*/

#if !defined(HASH_UTIL_KHUIEY3789Y98PDY89P3WAYU8W9PY)
#define HASH_UTIL_KHUIEY3789Y98PDY89P3WAYU8W9PY

#include "hash_map.h"

namespace ext
{

// returns the number of colliding pairs in the hash table
template <class HashType>
long hash_collisions(HashType const& h);

// returns the expectation value of the number of
// colliding pairs in the hash table.
// The number of pairs is (h.size() * (h.size() - 1)) / 2,
// and the probability of a pair colliding is 1 / h.bucket_count().
// Therefore the expectation value is (h.size() * (h.size() - 1)) / (2 * h.bucket_count()).
template <class HashType>
double hash_expected_collisions(HashType const& h);

// returns the efficiency of the hash function, defined as the ratio of
// the actual number of collisions to the expected number of collisions.
// Lower numbers are better.
template <class HashType>
inline
double hash_collision_ratio(HashType const& h)
{
   return hash_collisions(h) / hash_expected_collisions(h);
}

//
// implementation
//

template <class HashType>
long hash_collisions(HashType const& h)
{
   typedef typename HashType::size_type size_type;
   long Collide = 0;
   size_type Buckets = h.bucket_count();
   for (size_type b = 0; b < Buckets; ++b)
   {
      long BucketSize = h.bucket_size(b);
      Collide += (BucketSize * (BucketSize - 1)) / 2;
   }
   return Collide;
}

template <class HashType>
double hash_expected_collisions(HashType const& h)
{
   double Size = h.size();
   return Size * (Size - 1) / (h.bucket_count() * 2);
}

} // namespace ext

#endif
