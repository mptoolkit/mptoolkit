// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/randutil.h
//
// Copyright (C) 2016-2024 Ian McCulloch <ian@qusim.net>
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

// Simple functions for C++ random numbers.
// None of these functions are safe to call from multiple threads,
// except for crypto_rand(), which is thread-safe.
// The default seed is a constant -- applications should call one of the
// seed() functions before generating random numbers.
//
// If multiple independent streams of random numbers are required,
// then use the class random_stream.  It is safe to use different
// random_stream objects in different threads.
//
// 2024-02-08: Changed the default seeding to random (crypto) seed rather than deterministic

#if !defined(MPTOOLKIT_COMMON_RANDUTIL_H)
#define MPTOOLKIT_COMMON_RANDUTIL_H

#include <random>
#include <initializer_list>
#include <vector>

namespace randutil
{

// the random number generator.  This is the basic generator that is used
// to construct a specific distribution, which uses a 32-bit Mersenne Twister.
// It is an object that is declared here, but defined somewhere else
// (in the randutil.cpp implementation file) so we declare it as 'extern'
//
// We can also get an unsigned integer with
// unsigned x = u_rand();

extern std::mt19937 u_rand;

// A random generator for 'cryptographically secure' random numbers.
// These can be a limited resource, only use these sparingly, eg as
// a seed for a pseudo-random generator.  The seed() function uses this generator
// to seed the Mersenne Twister.

// return a single cryptographically secure random unsigned integer
unsigned crypto_rand();

// return an array of n cryptographically secure random unsigned integers
std::vector<unsigned> crypto_rand_vector(int n);

// functions to get random numbers from specific distributions

// a random integer in the closed interval [Min, Max]
int rand_int(int Min, int Max);

// returns a real number in the range [0,1)
double rand();

// returns a uniformly distributed real number, with mean 0, standard deviation 1
double randn();

//
// functions seed the generator
//

// re-seed the generator from 256 bits of 'cryptographically secure' random numbers
void seed();

// seed from a single integer
void seed(unsigned s);

// seed from an array
void seed(std::vector<unsigned> const& s);

// seed from a list
template <typename T>
void seed(std::initializer_list<T> s);

// returns the seed that was previously set
std::vector<unsigned> get_seed();

//
// Similar interface using a class
//

class random_stream
{
   public:
      // default constructor
      random_stream();

      // not copyable
      random_stream(random_stream const&) = delete;

      // not copy-assignable
      random_stream& operator=(random_stream const&) = delete;

      // movable
      random_stream(random_stream&&) = default;

      // move-assignable
      random_stream& operator=(random_stream&&) = default;

      // construct with an initial seed
      explicit random_stream(unsigned s);

      // construct with an initial seed from an array
      explicit random_stream(std::vector<unsigned> const& s);

      // construct with an initial seed from a list
      template <typename T>
      explicit random_stream(std::initializer_list<T> s);

      // the random number generator
      std::mt19937 u_rand;

      int rand_int(int Min, int Max);

      // returns a real number in the range [0,1)
      double rand();

      // returns a uniformly distributed real number, with mean 0, standard deviation 1
      double randn();

      // seed the generator from 256 bits of 'cryptographically secure' random numbers
      void seed();

      // seed from a single integer
      void seed(unsigned s);

      // seed from an array
      void seed(std::vector<unsigned> const& s);

      // seed from a list
      template <typename T>
      void seed(std::initializer_list<T> s);

      // returns the seed that was previously set
      std::vector<unsigned> get_seed();

   private:
      std::vector<unsigned> Seed;
      std::uniform_real_distribution<double> UniformDist;
      std::normal_distribution<double> NormalDist;
};

} // namespace

// include templates and definitions of inline functions
#include "randutil.cc"

#endif
