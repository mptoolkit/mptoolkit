// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/random.h
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

// Simple functions for C++ random numbers

#if !defined(MPTOOLKIT_COMMON_RANDUTIL_H)
#define MPTOOLKIT_COMMON_RANDUTIL_H

#include <random>

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
//
// We can also get a cryptographically random unsigned integer with
// unsigned x = crypto_rand();

extern std::random_device crypto_rand;

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

} // namespace

// include templates and definitions of inline functions
#include "randutil.cc"

#endif
