// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/randutil.cpp
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

#include "randutil.h"
#include <mutex>

namespace randutil
{

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

std::vector<unsigned> crypto_rand_vector(int n)
{
   std::lock_guard<std::mutex> guard(rd_mutex);
   std::vector<unsigned> Result;
   Result.reserve(n);
   for (int i = 0; i < n; ++i)
      Result.push_back(rd());
   return Result;
}

std::vector<unsigned> Seed = crypto_rand_vector(8);

// some awkwardness: the mt19937 requires an l-value reference, we can't initialize it from a temporary
std::seed_seq TempInitializer(Seed.begin(), Seed.end());

std::mt19937 u_rand(TempInitializer);

void seed()
{
   seed(crypto_rand_vector(8));
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
//   : Seed{1,2,3,4,5,6},//crypto_rand_vector(8)},
   : Seed{crypto_rand_vector(8)},
     UniformDist(0,1),
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
