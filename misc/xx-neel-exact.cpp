// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/xx-neel-exact.cpp
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

/*
  Calculates the exact solution of the time evolution of the XX model
  starting from an initial neel state.
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/bessel.hpp>

using namespace std;
using namespace boost::math;

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      cerr << "usage: xx-neel-exact <t>\n";
      return 1;
   }

   double t = boost::lexical_cast<double>(argv[1]);

   cout.precision(20);
   cout << setw(25) << t << ' ' << setw(25) << 0.5 * cyl_bessel_j(0, 2*t) << endl;
}
