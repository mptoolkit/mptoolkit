// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/statistics.h
//
// Copyright (C) 2001-2017 Ian McCulloch <ian@qusim.net>
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
  statistics.h

  A few statistics manipulations

  Created 2001-04-25 Ian McCulloch
  2012-02-04: Added size() function to moving_average
*/

#if !defined(STATISTICS_H_DSHR489UFHVFEIUNF48U398UCSAHU)
#define STATISTICS_H_DSHR489UFHVFEIUNF48U398UCSAHU

#include <iterator>
#include <deque>
#include <numeric>
#include <cstdlib>     // for abs(int)

namespace statistics
{

// distance is a useful stats function
using std::distance;
using std::size_t;

// Greatest Common Divisor using Euclid's algorithm for positive numbers
inline
int gcd_pos(int a, int b)
{
   if (a == b)
      return a;
   if (a > b)
      return gcd_pos(a-b, b);
   return gcd_pos(a, b-a);
}

inline
int gcd(int a, int b)
{
   return gcd_pos(std::abs(a),std::abs(b));
}

inline
int lcm(int a, int b)
{
   return std::abs(a*b) / gcd(a,b);
}

inline
int lcm(int a, int b, int c)
{
   return lcm(lcm(a,b),c);
}

inline
int lcm(int a, int b, int c, int d)
{
   return lcm(lcm(a,b),lcm(c,d));
}

// calculates the mean of a sequence.  Returns 0 if the sequence is empty
template <class FwdIter>
auto mean(FwdIter start, FwdIter end) -> decltype(*start * 1.0);

// calculates the variance with a given mean.  Formula is
// variance = sum_i (\bar{x} - x_i)^2 / (N - dof)
template <typename FwdIter, typename T>
double variance(FwdIter start, FwdIter end, T mean, int dof = 0);

// helper function to calculate the variance when the mean is not known.
// this is slower but numerically more stable than the one-pass formula
// <x^2> - <x>^2.
template <class FwdIter>
inline
double variance(FwdIter start, FwdIter end)
{
   return variance(start, end, mean(start, end), 1);
}

// calculates the chi2 value between two different sets,
// chi2 = sum (x_i - y_i)^2
template <class FwdIter1, class FwdIter2>
double chi2(FwdIter1 x_start, FwdIter1 x_end, FwdIter2 y_start);

// structure for data relating to a linear fit,
// of data points to the function y = m * x + c
struct linear_fit_result
{
   double m, c, chi2, variance, variance_m, variance_c;
   size_t N;
};

template <class FwdIter1, class FwdIter2, class FwdIter3>
linear_fit_result linear_fit(FwdIter1 x_start, FwdIter1 x_end,
                             FwdIter2 y_start,
                             FwdIter3 variance_start);

template <typename T>
class moving_average
{
   public:
      typedef T value_type;

      explicit moving_average(int Count) : Count_(Count) {}

      void push(T const& x)
      { Data.push_back(x); if (int(Data.size()) > Count_) Data.pop_front(); }

      void clear() { Data.clear(); }

      value_type value() const
      { return std::accumulate(Data.begin(), Data.end(), T()) / Data.size(); }

      value_type operator()() const
      { return this->value(); }

      // returns the actual number of items in the average
      std::size_t size() const
      { return Data->size(); }

   private:
      moving_average(); // not implemented

      int Count_;
      std::deque<value_type> Data;
};

template <typename T>
class moving_exponential
{
   public:
      typedef T value_type;

      moving_exponential() : Factor_(0.0) {}

      explicit moving_exponential(double RelaxationFactor) : Factor_(RelaxationFactor), Value(), Accum(0.0) {}

      void push(T const& x)
      { Value = (Value*Factor_) + x; Accum = (Accum*Factor_) + 1.0; }

      void clear() { Value = 0; Accum = 0; }

      value_type value() const
      { return (1.0 / Accum) * Value; }

      value_type operator()() const
      { return this->value(); }

   private:
      double Factor_;
      value_type Value;
      double Accum;
};

} // namespace

#include "statistics.cc"

#endif
