// -*- C++ -*- $Id$

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
