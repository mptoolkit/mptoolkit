// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/lanczos.h
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
//
// Simple Lanczos algorithm.
// No restarts, designed for applications where the number of iterations is small.
//
// Note on stopping criteria:
// Lanczos implementations are often fairly sloppy on how exactly to define the stopping criteria.
// The length of the residual vector is r = ||Ax - theta x}}
// This quantity is shift-invariant - ie, we can add a constant to A, and then theta
// will shift by the same constant, leaving r the same.
// However, using simply r < tol isn't a suitable criteria, because r isn't scale invariant.
// If we multiply A by some contant, then r also scales by the same constant.
//
// Often people use the criteria r/theta < tol
// which solves the scale-invariance problem but now the stopping criteria changes
// if a constant is added to the operator.  Bennani and Braconniery
// (Stopping criteria for eigensolvers, Maria Bennani and Thierry Braconniery,
// CERFACS Technical Report TR/PA/94/22) recommend using
// r/A_F < tol, where A_F is the Frobenius norm of A.
// On a shift by a constant k, the Frobenius norm changes by <~ sqrt(k), which only
// partially solves the problem.  Here we use the 'spectral diameter', being
// the difference between largest and smallest Ritz values (as approximating the
// true spectral diameter of the matrix).


#if !defined(LANCZOS_H_H2348975894UFP389P0)
#define LANCZOS_H_H2348975894UFP389P0

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>
#include <fstream>

double const LanczosBetaTol = 1E-14;

template <typename VectorType, typename MultiplyFunctor>
double Lanczos(VectorType& Guess, MultiplyFunctor MatVecMultiply, int& Iterations,
               double& Tol, int MinIter = 2, int Verbose = 0)
{
   std::vector<VectorType>     v;          // the Krylov vectors
   std::vector<VectorType>     Hv;         // the Krylov vectors

   LinearAlgebra::Matrix<double> SubH(Iterations+1, Iterations+1, 0.0);

   VectorType w = Guess;

   double Beta = norm_frob(w);
   CHECK(!std::isnan(Beta));
   // double OrigBeta = Beta;      // original norm, not used
   w *= 1.0 / Beta;
   v.push_back(w);

   w = MatVecMultiply(v[0]);
   Hv.push_back(w);
   SubH(0,0) = real(inner_prod(v[0], w));
   w -= SubH(0,0) * v[0];

   Beta = norm_frob(w);
   if (Beta < LanczosBetaTol)
   {
      if (Verbose > 0)
         std::cerr << "lanczos: immediate return, invariant subspace found, Beta="
                   << Beta << '\n';
      Guess = v[0];
      Iterations = 1;
      Tol = Beta;
      return SubH(0,0);
   }

   // It isn't meaningful to do a Krylov algorithm with only one matrix-vector multiply,
   // but we postpone the check until here to allow the corner case of
   // Iterations==1 but the algorithm converged in 1 step,
   // which might happen for example if the Hilbert space is 1-dimensional
   // (external checks that set Iterations=max(Iterations,Dimension) are not unreasonable).
   CHECK(Iterations > 1)("Number of iterations must be greater than one")(Iterations);

   for (int i = 1; i < Iterations; ++i)
   {
      SubH(i, i-1) = SubH(i-1, i) = Beta;
      w *= 1.0 / Beta;
      v.push_back(w);
      w = MatVecMultiply(v[i]);
      Hv.push_back(w);
      w -= Beta*v[i-1];
      SubH(i,i) = real(inner_prod(v[i], w));
      w -= SubH(i,i) * v[i];
      Beta = norm_frob(w);


      if (Beta < LanczosBetaTol)
      {
         // Early return, we can't improve over the previous energy and eigenvector
         if (Verbose > 0)
            std::cerr << "lanczos: early return, invariant subspace found, Beta="
                      << Beta << ", iterations=" << (i+1) << '\n';
         Iterations = i+1;
         LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i+1),
                                                LinearAlgebra::range(0,i+1));
         LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
         double Theta = EValues[0];    // smallest eigenvalue
         double SpectralDiameter = EValues[i] - EValues[0];  // largest - smallest
         VectorType y = M(0,0) * v[0];
         for (int j = 1; j <= i; ++j)
            y += M(0,j) * v[j];
         Tol = Beta / SpectralDiameter;
         Guess = y;
         return Theta;
      }

      // solution of the tridiagonal subproblem
      LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i+1),
                                             LinearAlgebra::range(0,i+1));
      if (std::isnan(M(0,0)))
      {
         std::ofstream Out("lanczos_debug.txt");
         Out << "NAN encountered in Lanczos\n"
             << "Beta=" << Beta << "\n\n"
             << "norm_frob(Guess)=" << norm_frob(Guess) << "\n\n"
             << "Guess=" << Guess << "\n\n"
             << "M=" << "\n\n"
             << "SubH=" << SubH << "\n\n";
         for (unsigned n = 0; n < v.size(); ++n)
         {
            Out << "V[" << n << "]=" << v[n] << "\n\n";
         }
      }

      LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
      double Theta = EValues[0];    // smallest eigenvalue
      double SpectralDiameter = EValues[i] - EValues[0];  // largest - smallest
      VectorType y = M(0,0) * v[0];
      for (int j = 1; j <= i; ++j)
         y += M(0,j) * v[j];

      // residual = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int j = 0; j < i; ++j)
         r += M(0,j) * Hv[j];

      double ResidNorm = norm_frob(r);

      if (ResidNorm < fabs(Tol * SpectralDiameter) && i+1 >= MinIter)
         //if (ResidNorm < Tol && i+1 >= MinIter)
      {
         if (Verbose > 0)
            std::cerr << "lanczos: early return, residual norm below threshold, ResidNorm="
                      << ResidNorm << ", SpectralDiameter=" << SpectralDiameter
                      << ", iterations=" << (i+1) << '\n';
         Iterations = i+1;
         Tol = ResidNorm / SpectralDiameter;
         Guess = y;
         return Theta;
      }

      if (i == Iterations-1) // finished?
      {
         Guess = y;
         Tol = -ResidNorm / SpectralDiameter;
         return Theta;
      }
   }

   PANIC("Should never get here");
   return -1.0;
}

// backwards-compatible version
template <typename VectorType, typename MultiplyFunctor>
double Lanczos(VectorType& Guess, MultiplyFunctor MatVecMultiply, int Iterations)
{
   int Iter = Iterations;
   double Tol = 1E-10;
   return Lanczos(Guess, MatVecMultiply, Iter, Tol);
}

#endif
