// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/randutil.cc
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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

#include <cmath>

namespace randutil
{

// implementation details
namespace detail
{
   extern std::uniform_real_distribution<double> UniformDist;
   extern std::normal_distribution<double> NormalDist;
} // namespace detail


inline
int rand_int(int Min, int Max)
{
   return int(std::floor(randutil::rand() * (Max-Min+1))) + Min;
}

// returns a real number in the range [0,1)
inline
double rand()
{
   return detail::UniformDist(u_rand);
}

// returns a uniformly distributed real number
inline
double randn()
{
   return detail::NormalDist(u_rand);
}

template <typename T>
void seed(std::initializer_list<T> s)
{
   seed(std::vector<unsigned>(s));
}

// class random_stream

inline
int random_stream::rand_int(int Min, int Max)
{
   return int(std::floor(this->rand() * (Max-Min+1))) + Min;
}

inline
double random_stream::rand()
{
   return UniformDist(u_rand);
}

inline
double random_stream::randn()
{
   return NormalDist(u_rand);
}

template <typename T>
void random_stream::seed(std::initializer_list<T> s)
{
   this->seed(std::vector<unsigned>(s));
}

} // namespace
