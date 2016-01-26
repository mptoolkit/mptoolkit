// -*- C++ -*- $Id$

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
