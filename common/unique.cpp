// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/hash.cpp
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

#include "unique.h"
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/utsname.h>
#include <algorithm>
#include <mutex>
#include <vector>

namespace ext
{

namespace
{

// 'salt' used for the get_unique function.
// On each call to get_unique(), GlobalSalt is used as part of the
// hash, and on exit GlobalSalt is set to the new hash value.
// GlobalSalt is protected by SaltMutex.
uint32_t GlobalSalt = 0xdeadbeef;

std::mutex SaltMutex;

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a
  structure that could supported 2x parallelism, like so:
      a -= b;
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c)                              \
{                                               \
  a -= b; a -= c; a ^= (c>>13);                 \
  b -= c; b -= a; b ^= (a<<8);                  \
  c -= a; c -= b; c ^= (b>>13);                 \
  a -= b; a -= c; a ^= (c>>12);                 \
  b -= c; b -= a; b ^= (a<<16);                 \
  c -= a; c -= b; c ^= (b>>5);                  \
  a -= b; a -= c; a ^= (c>>3);                  \
  b -= c; b -= a; b ^= (a<<10);                 \
  c -= a; c -= b; c ^= (b>>15);                 \
}

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

uint32_t hash_bytes(unsigned char const* k, size_t length, uint32_t initval = 0xdeadbeef)
{
   uint32_t a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;         /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
      b += (k[4] +((uint32_t)k[5]<<8) +((uint32_t)k[6]<<16) +((uint32_t)k[7]<<24));
      c += (k[8] +((uint32_t)k[9]<<8) +((uint32_t)k[10]<<16)+((uint32_t)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((uint32_t)k[10]<<24);
   case 10: c+=((uint32_t)k[9]<<16);
   case 9 : c+=((uint32_t)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((uint32_t)k[7]<<24);
   case 7 : b+=((uint32_t)k[6]<<16);
   case 6 : b+=((uint32_t)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((uint32_t)k[3]<<24);
   case 3 : a+=((uint32_t)k[2]<<16);
   case 2 : a+=((uint32_t)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}

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

uint32_t hash_buffer::hash() const
{
   return hash_bytes(&Buf[0], Buf.size());
}

uint32_t hash_buffer::hash(uint32_t Salt) const
{
   return hash_bytes(&Buf[0], Buf.size(), Salt);
}

void hash_buffer::push_back(void const* buf, size_t size)
{
   for (size_t i = 0; i < size; ++i)
   {
      Buf.push_back((static_cast<unsigned char const*>(buf))[i]);
   }
}

void hash_buffer::push_back(char n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(signed char n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(unsigned char n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(short n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(unsigned short n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(int n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(unsigned int n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(long n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(unsigned long n)
{
   push_back(&n, sizeof(n));
}

#if defined(USE_LONGLONG)
void hash_buffer::push_back(long long n)
{
   push_back(&n, sizeof(n));
}

void hash_buffer::push_back(unsigned long long n)
{
   push_back(&n, sizeof(n));
}
#endif // USE_LONGLONG

void hash_buffer::push_back(char const* str)
{
   push_back(str, strlen(str));
}

void hash_buffer::push_back(std::string const& str)
{
   push_back(str.c_str(), str.size());
}

} // namespace ext::unnamed

uint32_t get_unique()
{
   // form a hash of...
   hash_buffer HB;

   // the pid...
   HB.push_back(getpid());

   // the current time...
   struct timeval t;
   gettimeofday(&t, NULL);
   HB.push_back(t.tv_sec);
   HB.push_back(t.tv_usec);

   // the hostname...
   struct utsname u;
   uname(&u);
   HB.push_back(u.nodename);

   // and the 'salt'
   std::lock_guard<std::mutex> Lock(SaltMutex);
   GlobalSalt = HB.hash(GlobalSalt);

   return GlobalSalt;
}

} // namespace ext
