// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/randutil.cpp
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

#include "randutil.h"
#include <mutex>

namespace randutil
{

std::vector<unsigned> Seed{1,2,3,4,5,6};

// some awkwardness: the mt19937 requires an l-value reference, we can't initialize it from a temporary
std::seed_seq TempInitializer(Seed.begin(), Seed.end());

std::mt19937 u_rand(TempInitializer);

std::mutex rd_mutex;
std::random_device rd;

namespace detail
{
   std::uniform_real_distribution<double> UniformDist(0,1);
   std::normal_distribution<double> NormalDist;
} // namespace detail

unsigned crypto_rand()
{
   std::lock_guard<std::mutex> guard(rd_mutex);
   return rd();
}

void seed()
{
   seed({crypto_rand(), crypto_rand(), crypto_rand(), crypto_rand(), 
	 crypto_rand(), crypto_rand(), crypto_rand(), crypto_rand()});
}

// seed from a single integer
void seed(unsigned s)
{
   seed(std::vector<unsigned>(1,s));
}

// seed from an array
void seed(std::vector<unsigned> const& s)
{
   Seed = s;
   std::seed_seq SS(Seed.begin(), Seed.end());
   u_rand.seed(SS);
}

std::vector<unsigned> get_seed()
{
   return Seed;
}

// random_stream

random_stream::random_stream()
   : Seed{1,2,3,4,5,6},
     UniformDist(1,0),
     NormalDist()
{
   std::seed_seq TempInitializer(Seed.begin(), Seed.end());
   u_rand = std::mt19937(TempInitializer);
}

void
random_stream::seed()
{
   this->seed({crypto_rand(), crypto_rand(), crypto_rand(), crypto_rand(), 
	    crypto_rand(), crypto_rand(), crypto_rand(), crypto_rand()});
}

void
random_stream::seed(unsigned s)
{
   this->seed(std::vector<unsigned>(1,s));
}

void
random_stream::seed(std::vector<unsigned> const& s)
{
   Seed = s;
   std::seed_seq SS(Seed.begin(), Seed.end());
   u_rand.seed(SS);
}

} // namespace random
