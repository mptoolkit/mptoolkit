// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// utils/averagecorr.cpp
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
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include "common/math_const.h"
#include "common/stringutil.h"
#include <sstream>

using std::cout;
using std::setprecision;
using std::setw;
typedef std::complex<double> complex;

/*
  some theory:

  we calculate the k space correlation function by assuming that the averaged real space
  correlation function is reflection symmetric (ie an even function), which vanishes
  at +/- infinity.  Let F(x) be the real space correlation function.  F(x) for x >= 0
  is given by the program input, and we take F(x) for x < 0 to be F(-x) = F(x).
  Thus the discrete fourier transform is f(k) = sum_x F(x) * exp(i * pi * k * x / N).
  Here x,k ranges from -N+1 to +N.  The fourier transform of a real,even function 
  is itself real and even, so we only ever need to deal with x,k >= 0.  Thus the
  k space transform is f(0) = sum_x F(x), and f(x) = 2 sum_x cos(pi * k * x / N) for 0 < x < N.
*/

int main(int argc, char** argv)
{
   if (argc < 2 || argc > 3)
   {
      std::cerr << "usage: " << argv[0] << " <file> [<zero point value>]\n\n"
                << "expected file format:\n"
                << "    <x1> <x2> <real> <imag>\n\n"
                << "outputs:\n"
                << "    .x           real space\n"
                << "    .k           momentum space\n";
      return 1;
   }

   std::string filename = argv[1];

   typedef std::map<int, std::complex<double> > CorrDataType;

   CorrDataType CorrData;
   std::map<int, int> CorrCount;

   std::ifstream InFile(filename.c_str());
   int x, y;
   double CorrReal, CorrImag;
   int CorrLength = 0;
   while (InFile >> x >> y >> CorrReal >> CorrImag)
   {
      CorrData[std::abs(x-y)] += complex(CorrReal, CorrImag);
      ++CorrCount[std::abs(x-y)];
      CorrLength = std::max(abs(x-y), CorrLength);
   }

   std::vector<complex> RealSpace(CorrLength+1, 0.0);
   for (CorrDataType::const_iterator I = CorrData.begin(); I != CorrData.end(); ++I)
   {
      RealSpace[I->first] = I->second / double(CorrCount[I->first]);
   }

   // Calculate the Fourier Transform

   std::vector<complex> KSpace(CorrLength+1, 0.0);

   for (int k = 0; k < CorrLength+1; ++k)
   {
      KSpace[k] += RealSpace[0];  // * cos(0.0)
      for (int x = 1; x < CorrLength+1; ++x)
      {
         KSpace[k] += 2.0 * RealSpace[x] * cos(double(math_const::pi * k * x) / (CorrLength+1));
      }
   }

   std::ofstream RealOut((filename + ".x").c_str());
   std::ofstream KOut((filename + ".k").c_str());

   RealOut << setprecision(12);
   KOut << setprecision(12);

   // write only the points of the real space correlation that actually exist
   for (CorrDataType::const_iterator I = CorrData.begin(); I != CorrData.end(); ++I)
   {
      RealOut << setw(3) << I->first << ' ' << RealSpace[I->first].real() << ' ' << RealSpace[I->first].imag() << '\n';
   }

   // write the k-space correlation
   for (int i = 0; i < CorrLength+1; ++i)
   {
      KOut << setw(14) << (double(i) / CorrLength) << ' ' << KSpace[i].real() << ' ' << KSpace[i].imag() << '\n';
   }

   return 0;
}

   
