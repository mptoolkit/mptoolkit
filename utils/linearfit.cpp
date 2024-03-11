// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// utils/linearfit.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include <iostream>
#include "common/statistics.h"
#include <math.h>
#include <vector>

using namespace std;
using namespace statistics;

int main()
{
   double x, y, dy;
   vector<double> X, Y, DY;
   cout.precision(20);
   while (cin >> x >> y >> dy)
   {
      X.push_back(x);
      Y.push_back(y);
      DY.push_back(dy);  // DY is the variance, not the standard deviation
   }

   linear_fit_result R = linear_fit(X.begin(), X.end(), Y.begin(), DY.begin());

   cout << R.c << '\n' << ' ' << sqrt(R.variance_c) << '\n';
}
