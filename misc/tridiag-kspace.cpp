// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/tridiag-kspace.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
#include <cmath>
#include "common/math_const.h"
#include <boost/lexical_cast.hpp>
#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"

using namespace LinearAlgebra;

int main(int argc, char** argv)
{
   if (argc != 1)
   {
      std::cerr << "usage: tridiag-kspace < file-of-hopping-energy-pairs\n";
      return 1;
   }

   std::vector<double> Hop, Epsilon;
   Epsilon.push_back(0.0); // the energy of the center of the star, take it to be 0
   {
      double h,e;
      while (std::cin >> e >> h)
      {
         Hop.push_back(h);
         Epsilon.push_back(e);
      }
   }

   int const d = Epsilon.size();
   Matrix<double> M(d, d, 0.0);
   for (unsigned i = 0; i < Hop.size(); ++i)
   {
      M(i,i) = Epsilon[i];
      M(0,i+1) = M(i+1,0) = Hop[i];
   }
   M(d-1,d-1) = Epsilon[d-1];

   TRACE(M);

   Matrix<std::complex<double> > Mc = M;
   Matrix<double> T = real(TridiagonalizeHermitian(Mc));

   TRACE(T);

   for (unsigned i = 0; i < Hop.size(); ++i)
   {
      Epsilon[i] = T(i,i);
      Hop[i] = T(i,i+1);
   }
    Epsilon[d-1] = T(d-1,d-1);

   std::cout << "Energies:\n" << Epsilon
             << "\n\nHoppings:\n" << Hop << '\n';

}
