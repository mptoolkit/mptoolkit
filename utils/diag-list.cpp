// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// utils/diag-list.cpp
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

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include <complex>

using namespace LinearAlgebra;

int main()
{
   std::vector<std::complex<double> > v;
   std::complex<double> z;
   while (std::cin >> z)
      v.push_back(z);

   int n = int(std::sqrt(v.size()));
   Matrix<std::complex<double> > M(n, n);
   std::vector<std::complex<double> >::const_iterator vI = v.begin();
   for (iterator<Matrix<std::complex<double> > >::type I = iterate(M); I; ++I)
      for (inner_iterator<Matrix<std::complex<double> > >::type J = iterate(I); J; ++J)
         *J = *vI++;

   std::cout.precision(12);
   std::cout << "Size: " << v.size() << '\n';
   std::cout << "Eigenvalues:\n";
   std::cout << EigenvaluesHermitian(M) << '\n';
   DiagonalizeHermitian(M);
   std::cout << "Eigenvectors:\n";
   std::cout << M << std::endl;
}
