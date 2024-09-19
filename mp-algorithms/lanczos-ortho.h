// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/lanczos-ortho.h
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
//
//  Version of lanczos.h that forces the Krylov subspace to be orthogonal to a set of
//  user-supplied vectors.
//

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

template <typename VectorType, typename MultiplyFunctor>
double Lanczos(VectorType& Guess, MultiplyFunctor MatVecMultiply, int& Iterations,
               std::vector<VectorType> const& Ortho, bool UseDGKS = false)
{
   std::vector<VectorType>     v(Iterations+1);         // the Krylov vectors
   std::vector<double>         Alpha(Iterations+1);
   std::vector<double>         Beta(Iterations+1);

   VectorType r = Guess;
   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      r -= inner_prod(Ortho[j], r) * Ortho[j];
   }
   if (UseDGKS)
   {
      for (std::size_t j = 0; j < Ortho.size(); ++j)
      {
         r -= inner_prod(Ortho[j], r) * Ortho[j];
      }
   }

   Beta[0] = norm_frob(r);

   double tol = 1E-11;

   int i;
   for (i = 1; i <= Iterations; ++i)
   {
      v[i] = (1.0 / Beta[i-1]) * r;
      r = MatVecMultiply(v[i]);
      for (std::size_t j = 0; j < Ortho.size(); ++j)
      {
         r -= inner_prod(Ortho[j], r) * Ortho[j];
      }
      // DGKS-like correction
      if (UseDGKS)
      {
         for (std::size_t j = 0; j < Ortho.size(); ++j)
         {
            r -= inner_prod(Ortho[j], r) * Ortho[j];
         }
      }
      if (i > 1)
         r = r - Beta[i-1] * v[i-1];
      Alpha[i] = real(inner_prod(v[i], r));
      r = r - Alpha[i] * v[i];
      Beta[i] = norm_frob(r);

      if (Beta[i] <= tol)
      {
        ++i;
        break;
      }
   }
   Iterations = i-1;

   LinearAlgebra::Matrix<double> M(Iterations, Iterations, 0.0);
   M(0,0) = Alpha[1];
   for (int k = 1; k < Iterations; ++k)
   {
      M(k,k) = Alpha[k+1];
      M(k,k-1) = M(k-1,k) = Beta[k];
   }

   LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);

   DEBUG_TRACE(EValues);

   // calculate the ground-state Lanczos vector
   Guess = M(0,0) * v[1];
   for (int k = 1; k < i-1; ++k)
   {
      Guess += M(0,k) * v[k+1];
   }

   // Normalize the result vector
   Guess *= Beta[0] / norm_frob(Guess);
   DEBUG_TRACE(norm_frob(Guess));

   //   TRACE(norm_frob_sq(Guess));
   //   TRACE(Guess);

   return EValues[0];
}
