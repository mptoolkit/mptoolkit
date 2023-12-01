// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testhaar.cpp
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

#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/eigen.h"

using namespace LinearAlgebra;

// A random unitary matrix, distributed according to the Haar measure
inline
Matrix<std::complex<double> > random_unitary(size_type Size)
{
   Matrix<std::complex<double> > M = nrandom_matrix<std::complex<double> >(Size, Size);
   Matrix<std::complex<double> > R = M;
   Matrix<std::complex<double> > Q = QR_FactorizeFull(R);
   return Q;
   Matrix<std::complex<double> > Lambda(Size, Size, 0.0);
   for (int i = 0; i < Size; ++i)
   {
      Lambda(i,i) = R(i,i) / norm_frob(R(i,i));
   }
   return Q*Lambda;
}

Vector<double> GetEigenvalues()
{
   Matrix<std::complex<double> > x = random_unitary(64);
   Matrix<std::complex<double> > psi(2,32);
   for (int i = 0; i < 32; ++i)
   {
      psi(0,i) = x(0,i);
      psi(1,i) = x(0,i+32);
   }
   Matrix<std::complex<double> > rho = psi * herm(psi);
   return EigenvaluesHermitian(rho);
}

int main()
{
   int NumBuckets = 100;
   int NumRep = 5000;
   std::vector<int> Hist(NumBuckets, 0);

   for (int i = 0; i < NumRep; ++i)
   {
      Vector<double> e = GetEigenvalues();
      ++Hist[int(e[0]*NumBuckets)];
      ++Hist[int(e[1]*NumBuckets)];
   }

   for (int i = 0; i < NumBuckets; ++i)
   {
      std::cout << std::setw(20) << ((i+0.5)/NumBuckets) << ' ' << Hist[i] << '\n';
   }
}
