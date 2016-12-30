// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/spin1chain-minenergy.cpp
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
#include <iostream>
#include <cmath>
#include "common/math_const.h"

using namespace std;

int main()
{
   int const L = 480;

   double const ThetaDenom = 400.0 / (2*math_const::pi);
   int const ThetaMin = 0;
   int const ThetaMax = 499;

   double const BoundaryFalloff = 0.7;
   double BilinearEnergy[3] = {-2, -1, 1};
   double BiquadEnergy[3] = {4,1,1};

   double BilinearBoundaryEnergy[2] = {-1,0.5};
   double BiquadBoundaryEnergy[2] = {1,0.25};

   cout.precision(14);
   for (int i = ThetaMin; i <= ThetaMax; ++i)
   {
      double BondEnergy[3];
      for (int j = 0; j < 3; ++j)
      {
         BondEnergy[j] = cos(double(i)/ThetaDenom) * BilinearEnergy[j]
            + sin(double(i)/ThetaDenom) * BiquadEnergy[j];
      }
      double BondMin = std::min(BondEnergy[0], std::min(BondEnergy[1], BondEnergy[2]));

      double BoundaryEnergy[2];
      for (int j = 0; j < 2; ++j)
      {
         BoundaryEnergy[j] = cos(double(i)/ThetaDenom) * BilinearBoundaryEnergy[j]
            + sin(double(i)/ThetaDenom) * BiquadBoundaryEnergy[j];
      }
      double BoundaryMin = std::min(BoundaryEnergy[0], BoundaryEnergy[1]);

      double E = BoundaryMin * BoundaryFalloff * 2 + (L-2)*BondMin;
      std::cout << i << ' ' << E << '\n';
   }
}
