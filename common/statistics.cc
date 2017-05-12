// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/statistics.cc
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

namespace statistics
{

template <class FwdIter>
auto mean(FwdIter start, FwdIter end) -> decltype(*start * 1.0)
{
   if (start == end)
      return decltype(*start * 1.0){};

   size_t n = 1;
   decltype(*start * 1.0) x(*start);
   ++start;
   while (start != end)
   {
      x += *start;
      ++n;
      ++start;
   }
   return x / double(n);
}

template <typename FwdIter, typename T>
double variance(FwdIter start, FwdIter end, T mean, int dof)
{
   using std::norm;
   size_t n = 0;
   double v = 0;
   while (start != end)
   {
      v += norm(mean - *start);
      ++n;
      ++start;
   }
   return v / (n - dof);
}

template <class FwdIter1, class FwdIter2, class FwdIter3>
linear_fit_result linear_fit(FwdIter1 x_start, FwdIter1 x_end, FwdIter2 y_start, FwdIter3 variance_start)
{
   linear_fit_result result;
   result.N = 0;
   double delta = 0;
   double X = 0, Y = 0, X2 = 0, Y2 = 0, XY = 0, E = 0;
   // X = sum_i {x_i / var_i}
   // Y = sum_i {y_i / var_i}
   // X2 = sum_i {x_i^2 / var_i}
   // Y2 = sum_i {y_i^2 / var_i}
   // XY = sum_i {x_i y_i / var_i}
   // E = sum_i {1 / var_i}
   FwdIter1 x = x_start;
   FwdIter2 y = y_start;
   FwdIter3 v = variance_start;
   while (x != x_end)
   {
      double e = 1.0 / (*v);  // e = 1.0/variance
      X += (*x) * e;
      Y += (*y) * e;
      X2 += (*x) * (*x) * e;
      Y2 += (*y) * (*y) * e;
      XY += (*x) * (*y) * e;
      E += e;
      ++result.N;
      ++x;
      ++y;
      ++v;
   }

   delta = E * X2 - X * X;
   result.m = (1 / delta) * (E * XY - X * Y);
   result.c = (1 / delta) * (X2 * Y - X * XY);
   result.variance_m = (1 / delta) * E;
   result.variance_c = (1 / delta) * X2;

   // calculate chi2 and the variance
   x = x_start;
   y = y_start;
   v = variance_start;
   result.chi2 = 0;
   result.variance = 0;
   while (x != x_end)
   {
      double t = (*y - (result.m * (*x) + result.c));
      result.chi2 += t * t / (*v);
      result.variance += t * t;
      ++x;
      ++y;
      ++v;
   }
   result.variance /= (result.N - 2);  // 2 dof in the calculation, m & c

   return result;
}

template <class FwdIter1, class FwdIter2>
double chi2(FwdIter1 x_start, FwdIter1 x_end, FwdIter2 y_start)
{
   double result = 0;
   while (x_start != x_end)
   {
      double t = (*x_start - *y_start);
      result += t * t;
      ++x_start;
      ++y_start;
   }
   return result;
}

} // namespace
