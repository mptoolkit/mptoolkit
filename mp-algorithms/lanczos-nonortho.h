// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/lanczos-nonortho.h
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

// lanczos function for a non-orthogonal basis.  This means that we
// additionally supply functors for the inner_prod and the norm_frob.

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include "lanczos.h"
#include <iostream>
#include <cmath>

#if !defined(LANCZOS_NONORTHO_H_FDJIH834957UY89JP89)
#define LANCZOS_NONORTHO_H_FDJIH834957UY89JP89

template <typename VectorType, typename MultiplyFunctor,
          typename InnerProdType, typename NormFrobType>
double Lanczos(VectorType& Guess, MultiplyFunctor MatVecMultiply,
               InnerProdType InnerProd,
               NormFrobType NormFrob, int& Iterations,
               double& Tol, bool Verbose = false)
{
   std::vector<VectorType>     v;          // the Krylov vectors
   std::vector<VectorType>     Hv;         // the Krylov vectors

   LinearAlgebra::Matrix<double> SubH(Iterations+1, Iterations+1, 0.0);

   VectorType w = Guess;

   double Beta = NormFrob(w);
   // double OrigBeta = Beta;      // original norm, not used
   w *= 1.0 / Beta;
   v.push_back(w);

   w = MatVecMultiply(v[0]);
   Hv.push_back(w);
   SubH(0,0) = real(InnerProd(v[0], w));
   w -= SubH(0,0) * v[0];

   Beta = NormFrob(w);
   if (Beta < LanczosBetaTol)
   {
      if (Verbose)
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
      SubH(i,i) = real(InnerProd(v[i], w));
      w -= SubH(i,i) * v[i];
      Beta = NormFrob(w);

      // solution of the tridiagonal subproblem
      LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i+1),
                                             LinearAlgebra::range(0,i+1));
      LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
      double Theta = EValues[0];    // smallest eigenvalue
      VectorType y = M(0,0) * v[0];
      for (int j = 1; j <= i; ++j)
         y += M(0,j) * v[j];

      // residual = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int j = 0; j < i; ++j)
         r += M(0,j) * Hv[j];

      double ResidNorm = NormFrob(r);

      if (ResidNorm < fabs(Tol * Theta))
      {
         if (Verbose)
            std::cerr << "lanczos: early return, residual norm below threshold, ResidNorm="
                      << ResidNorm << ", iterations=" << (i+1) << '\n';
         Iterations = i+1;
         Tol = ResidNorm/fabs(Theta);
         Guess = y;
         return Theta;
      }

      if (Beta < LanczosBetaTol)
      {
         if (Verbose)
            std::cerr << "lanczos: early return, invariant subspace found, Beta="
                      << Beta << ", iterations=" << (i+1) << '\n';
         Iterations = i+1;
         Tol = Beta;
         Guess = y;
         return Theta;
      }

      if (i == Iterations-1) // finished?
      {
         Guess = y;
         Tol = -ResidNorm/fabs(Theta);
         return Theta;
      }
   }

   PANIC("Should never get here");
   return -1.0;
}

#endif
