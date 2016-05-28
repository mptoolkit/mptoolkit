// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/continued-fraction.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <boost/lexical_cast.hpp>

using namespace std;

complex<double>
evaluate_fraction(complex<double> const& z, vector<complex<double> > const& c)
{
   complex<double> x = 0;
   for (int i = c.size()-1; i >= 0; --i)
   {
      x = c[i] / (z + x);
   }
   return x;
}

int main(int argc, char** argv)
{
   if (argc == 1)
   {
      cerr << "usage: continued-fraction x [ x ... ] < list_of_coefficients_cn\n"
           << "\nEvaluates the continued fraction c(z) = c0 / (z + c1 / (z + c2 / (z + .... )))\n";
      return 1;
   }

   vector<complex<double> > c;
   double x;
   while (cin >> x)
      c.push_back(x);

   cout.precision(16);

   for (int i = 1; i < argc; ++i)
   {
      complex<double> z = boost::lexical_cast<std::complex<double> >(argv[i]);
      complex<double> x = evaluate_fraction(z, c);
      std::cout << std::setw(20) << z << "   " << std::setw(20) << x << '\n';
   }
}
