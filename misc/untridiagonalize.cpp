// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// misc/untridiagonalize.cpp
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
      std::cerr << "usage: untridiagonalize < bath-file\n";
      return 1;
   }

   std::vector<double> Hop, Epsilon;
   {
      double h, e;
      // read the first epsilon term, after that we get (epsilon,hop) pairs
      std::cin >> e;
      Epsilon.push_back(e);
      while (std::cin >> h >> e)
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
      M(i,i+1) = M(i+1,i) = Hop[i];
   }
   M(d-1,d-1) = Epsilon[d-1];

   Matrix<double> Mtail = M(range(1,d), range(1,d));
   Matrix<double> U = Mtail;
   Vector<double> Eigen = DiagonalizeHermitian(U);
   
   Matrix<double> Ufull(d,d,0.0);
   Ufull(0,0) = 1;
   Ufull(range(1,d), range(1,d)) = U;

   Matrix<double> Result = Ufull * M * trans(Ufull);

   std::cout << "Energies:\n" << direct_sum(M(0,0), Eigen)
             << "\n\nHoppings:\n" << Result(0, all)[range(1,size2(Result))] << '\n';

}
