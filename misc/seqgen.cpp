// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// misc/seqgen.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
#include "linearalgebra/scalarmatrix.h"

using namespace LinearAlgebra;

inline
Matrix<std::complex<double> > random_unitary(size_type Size)
{
   Matrix<std::complex<double> > M = nrandom_matrix<std::complex<double> >(Size, Size);
   Matrix<std::complex<double> > R = M;
   Matrix<std::complex<double> > Q = QR_Factorize(R);
   Matrix<std::complex<double> > Lambda(Size, Size, 0.0);
   for (int i = 0; i < Size; ++i)
   {
      Lambda(i,i) = R(i,i) / norm_frob(R(i,i));
   }
   return Q*Lambda;
}

inline
Matrix<std::complex<double> > random_unitary(size_type Size1, size_type Size2)
{
   int sz = std::max(Size1, Size2);
   Matrix<std::complex<double> > Mat = random_unitary(sz);
   return Mat(range(0,Size1), range(0,Size2));
}

typedef Matrix<std::complex<double> > ComplexMatrix;

Vector<double> GetEigen(int const D, int const Chi, int const L)
{
   ComplexMatrix U = random_unitary(D*Chi, Chi);

   Vector<ComplexMatrix> A(2);
   A[0] = U(range(0,Chi), all);
   A[1] = U(range(Chi,2*Chi), all);

   ComplexMatrix RightBoundary = random_unitary(Chi,1);
   ComplexMatrix LeftBoundary = random_unitary(1,Chi);

   ComplexMatrix X = RightBoundary * herm(RightBoundary);
   for (int i = 0; i < L-1; ++i)
   {
      ComplexMatrix Xp(Chi, Chi, 0.0) ;
      for (int s = 0; s < D; ++s)
      {
         Xp += A[s] * X * herm(A[s]);
      }
      X = Xp;
   }

   // make the DM
   ComplexMatrix rho(D,D, 0.0);
   for (int s = 0; s < D; ++s)
   {
      for (int t = 0; t < D; ++t)
      {
         rho(s,t) = ComplexMatrix(LeftBoundary * A[s] * X * herm(A[t]) * herm(LeftBoundary))(0,0);
      }
   }

   TRACE(trace(rho));
   rho *= 1.0 / trace(rho);

   return EigenvaluesHermitian(rho);
}

void check(int const D, int const Chi, int const L)
{
   ComplexMatrix U = random_unitary(D*Chi, Chi);

   Vector<ComplexMatrix> A(2);
   A[0] = U(range(0,Chi), all);
   A[1] = U(range(Chi,2*Chi), all);

   ComplexMatrix RightBoundary = random_unitary(Chi,1);
   ComplexMatrix LeftBoundary = random_unitary(1,Chi);

   ComplexMatrix X = ScalarMatrix<double>(Chi, 1.0);

   for (int i = 0; i < L-1; ++i)
   {
      ComplexMatrix Xp(Chi, Chi, 0.0) ;
      for (int s = 0; s < D; ++s)
      {
         Xp += herm(A[s]) * X * A[s];
      }
      X = Xp;
   }
   TRACE(X);
}

int main()
{
   int const D = 2;
   int const Chi = 3;
   int const L = 4;

   check(D, Chi, L);
   PANIC("here");

   int NumBuckets = 100;
   int NumRep = 5000;
   std::vector<int> Hist(NumBuckets, 0);

   for (int i = 0; i < NumRep; ++i)
   {
      std::cout << i << '\n';
      Vector<double> e = GetEigen(D, Chi, L);
      ++Hist[int(e[0]*NumBuckets)];
      ++Hist[int(e[1]*NumBuckets)];
   }

   for (int i = 0; i < NumBuckets; ++i)
   {
      std::cout << std::setw(20) << ((i+0.5)/NumBuckets) << ' ' << Hist[i] << '\n';
   }
}
