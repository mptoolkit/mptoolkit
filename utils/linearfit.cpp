// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// utils/linearfit.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
