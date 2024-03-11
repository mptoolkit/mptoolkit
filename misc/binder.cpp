// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/binder.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
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

#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <cmath>

using namespace std;

int main(int argc, char** argv)
{
   if (argc != 7)
   {
      cerr << "usage: binder <1> <2> <3> <4> <corr> <scale>\n";
      return 1;
   }

   // the 4th order moments
   double t4_1 = boost::lexical_cast<double>(argv[1]);
   double t4_2 = boost::lexical_cast<double>(argv[2]);
   double t4_3 = boost::lexical_cast<double>(argv[3]);
   double t4_4 = boost::lexical_cast<double>(argv[4]);
   double xi = boost::lexical_cast<double>(argv[5]);
   double Scale = boost::lexical_cast<double>(argv[6]);

   // convert the components of the moments into the cumulants
   double k4 = t4_1;
   double k1 = sqrt(sqrt(t4_4));
   double k2 = t4_3 / (6.0 * sqrt(t4_4));
   double k3 = (t4_2 - 3.0*k2*k2)/(4.0*k1);

   // the 2nd order moments
   double t2_1 = k2;
   double t2_2 = sqrt(t4_4);

   // Effective length
   double L = Scale * xi;

   double m4 = (((t4_4*L + t4_3)*L + t4_2)*L + t4_1)*L;
   double m2 = (t2_2*L + t2_1)*L;

   // the Binder cumulant
   double U = 1.0 - m4 / (3.0 * m2*m2);

   cout.precision(14);
   cout << U << endl;
}
