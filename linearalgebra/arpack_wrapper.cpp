// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/arpack_wrapper.cpp
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

#include "arpack_wrapper.h"
#include <vector>
#include <iostream>

namespace LinearAlgebra
{

typedef std::complex<double> complex;

void swap_vectors(int n, complex* x, complex* y)
{
   std::vector<complex> Temp(x, x+n);
   fast_copy(y, y+n, x);
   fast_copy(&Temp[0], &Temp[0]+n, y);
}

void
MatchEigenvectors(int n,
                  Vector<std::complex<double> >& LeftValues,
                  std::vector<std::complex<double> >& LeftVectors,
                  Vector<std::complex<double> >& RightValues,
                  std::vector<std::complex<double> >& RightVectors, double tol)
{
   CHECK_EQUAL(LeftValues.size(), RightValues.size());

   for (unsigned i = 0; i < LeftValues.size(); ++i)
   {
      // find the right eigenvalue closest to LeftValues[i]
      unsigned Pos = i;
      double Dist = norm_frob(LeftValues[i] - RightValues[Pos]);
      for (unsigned j = i+1; j < RightValues.size(); ++j)
      {
         double jDist = norm_frob(LeftValues[i] - RightValues[j]);
         if (jDist < Dist)
         {
            Pos = j;
            Dist = jDist;
         }
      }
      if (Dist > tol)
         std::cerr << "MatchEigenvalues: warning: left & right eigenvalues differ by > tol.  Left="
                   << LeftValues[i] << ", Right=" << RightValues[Pos] << ", delta=" << Dist << '\n';
      std::swap(RightValues[i], RightValues[Pos]);
      swap_vectors(n, &RightVectors[n*i], &RightVectors[n*Pos]);
   }
}

} // namespace LinearAlgebra
