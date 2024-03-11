// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/hash.h
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
  hash.h

  Portable hashing.

  Created 2002-11-17 Ian McCulloch
*/

#if !defined(HASH_H_HDSFJKHUIRY879YHF8HQ38HF8OH389UWEO)
#define HASH_H_HDSFJKHUIRY879YHF8HQ38HF8OH389UWEO

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif
#include "common/inttype.h"
#include <string>
#include <vector>
#include <string.h>

namespace ext
{

using inttype::uint32_t;

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

uint32_t hash_bytes(unsigned char const* k, size_t length, uint32_t initval = 0xdeadbeef);

/*
  My comment (Ian McCulloch): By default, hash table sizes are primes, and expansion
  is by factor ~ 2.  The list of table sizes is given by FindNextPrime() which has a table of
  prime numbers taken from SGI's STL library.  If you define HASHTABLE_POWER2, the hash table
  size is instead always a power 2, and the modulo operation is replaced by a bit-wise 'and'
  with a mask.

  In benchmark tests on a pentium3, defining HASHTABLE_POWER2 gives a measurable but
  very small speedup - negligable for most practical purposes.  For general use, it is
  probably better to default to NOT defining HASHTABLE_POWER2, which gives better
  behaviour with poor hash functions.

  On a pentium3, the hash_map implementation benchmarks as slightly slower than
  the hash_map extension in the GNU libstdc++ library (which is itself derived
  from the SGI version).  The difference is presumably in
  the list handling, as the GNU version uses a hand-rolled singly-linked list,
  versus a std::list.  For the particular benchmark used (find() performance with
  a hash_map<int, int>), the speed difference is even smaller than the effect of
  HASHTABLE_POWER2, i.e. with HASHTABLE_POWER2 this version becomes microscopically
  faster than GNU.  The performance is much more erratic than GNU though, I don't know why.
*/

// default hash functor.  Specializations are defined for builtin integral types,
// char*, char const* (assumed to be C-style strings), std::string,
// void* and void const* (hashes on the value of the pointer).

template <typename Key> struct hash;

// class for calculating a hash from arbitary data.
// Possibly this should be templated on a hash function?
// Currently we use the hash_bytes() function.
struct hash_buffer
{
   public:
      hash_buffer() {}

      void push_back(char n);
      void push_back(signed char n);
      void push_back(unsigned char n);
      void push_back(short n);
      void push_back(unsigned short n);
      void push_back(int n);
      void push_back(unsigned int n);
      void push_back(long n);
      void push_back(unsigned long n);
#if defined(USE_LONGLONG)
      void push_back(long long n);
      void push_back(unsigned long long n);
#endif

      void push_back(char const* str);
      void push_back(std::string const& str);

      void push_back(void const* buf, size_t size);

      uint32_t hash() const;

      uint32_t hash(uint32_t Salt) const;

   private:
      std::vector<unsigned char> Buf;
};

// hash specializations

template <>
struct hash<char>
{
   size_t operator()(char i) const { return size_t(i); }
};

template <>
struct hash<unsigned char>
{
   size_t operator()(unsigned char i) const { return size_t(i); }
};

template <>
struct hash<signed char>
{
   size_t operator()(signed char i) const { return size_t(i); }
};

template <>
struct hash<int>
{
   size_t operator()(int i) const { return size_t(i); }
};

template <>
struct hash<unsigned int>
{
   size_t operator()(unsigned int i) const { return size_t(i); }
};

template <>
struct hash<short>
{
   size_t operator()(short i) const { return size_t(i); }
};
template <>
struct hash<unsigned short>
{
   size_t operator()(unsigned short i) const { return size_t(i); }
};

// this should do something different if sizeof(long) > sizeof(size_t) ???
template <>
struct hash<long>
{
   size_t operator()(long i) const { return size_t(i); }
};
template <>
struct hash<unsigned long>
{
   size_t operator()(unsigned long i) const { return size_t(i); }
};

#if defined(USE_LONGLONG)
template <>
struct hash<long long>
{
   size_t operator()(long long i) const { return size_t(i); }
};

template <>
struct hash<unsigned long long>
{
   size_t operator()(unsigned long long i) const { return size_t(i); }
};
#endif

template <>
struct hash<char const*>
{
  size_t operator()(char const* str) const
  {
    return hash<size_t>()(hash_bytes(reinterpret_cast<unsigned char const*>(str), strlen(str)));
  }
};

template <>
struct hash<char*>
{
  size_t operator()(char* str) const
  {
    return hash<size_t>()(hash_bytes(reinterpret_cast<unsigned char*>(str), strlen(str)));
  }
};

template <>
struct hash<std::string>
{
  size_t operator()(std::string const& str) const
  {
    return hash<size_t>()(hash_bytes(reinterpret_cast<unsigned char const*>(str.data()), str.size()));
  }
};

template <>
struct hash<void const*>
{
   size_t operator()(void const* p) const { return reinterpret_cast<size_t>(p); }
};

template <>
struct hash<void*>
{
   size_t operator()(void* p) const { return reinterpret_cast<size_t>(p); }
};

// Returns a hash constructed such that the return value is unique with a probability
// of 1 in 2^32.  The data hashed includes the current time from gettimeofday(),
// the process id, the hostname, and the value from the previous call.
uint32_t get_unique();

// The default version of get_unique() uses a mutex on the value of the previous call.
// This means that it can't be used before dynamic initialization.  The
// no_mutex version doesn't grab the mutex, but is not threadsafe.
uint32_t get_unique_no_mutex();

namespace Private
{
// this is some stuff used by the hashed containers.  The implementation goes in
// hash.cpp, since otherwise we would need another .cpp file especially for this.
// Thus, the definition goes here.

unsigned long FindNextPrime(unsigned long p);

double const DefaultMaxLoadFactor = 0.8;

extern size_t const MinBuckets;  // minimum number of buckets, equal to the smallest prime in FindNextPrime()

#if defined(HASHTABLE_POWER2)
extern size_t const MinMask;
#endif

} // namespace ext::Private

} // namespace ext

#endif
