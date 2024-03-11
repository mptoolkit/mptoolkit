// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/lanczos-exponential.h
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

#include "linearalgebra/eigen.h"
#include "linearalgebra/exponential.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;

template <typename VectorType, typename MultiplyFunctor>
VectorType LanczosExponential(VectorType const& x,
                              MultiplyFunctor MatVecMultiply, int& Iterations,
                              std::complex<double> const& Theta,
                              double& ETol)
{
   typedef std::complex<double> complex;
   std::vector<VectorType>     v(Iterations+1);         // the Krylov vectors
   std::vector<double>         Alpha(Iterations+1);
   std::vector<double>         Beta(Iterations+1);

   VectorType r = x;
   Beta[0] = norm_frob(r);

   double const BetaTol = 1E-14;
   double Error = 1;

   LinearAlgebra::Matrix<double> M(Iterations+1, Iterations+1, 0.0);

   LinearAlgebra::Vector<complex> Coefficients(1, 1.0);

   int i;
   for (i = 1; i <= Iterations && Error > ETol; ++i)
   {
      v[i] = (1.0 / Beta[i-1]) * r;
      r = MatVecMultiply(v[i]);
      r = r - Beta[i-1] * v[i-1];
      Alpha[i] = real(inner_prod(v[i], r));
      r = r - Alpha[i] * v[i];
      Beta[i] = norm_frob(r);

      M(i-1,i-1) = Alpha[i];
      M(i,i-1) = M(i-1,i) = Beta[i];

      if (i == 1)
         Error = 1.0;
      else
      {
         LinearAlgebra::Matrix<complex> Msub = Theta * M(range(0,i), range(0,i));
         LinearAlgebra::Matrix<complex> X = LinearAlgebra::Exponentiate(1.0, Msub);
         //Error =  LinearAlgebra::norm_2(X(i-1,0));
         //DEBUG_TRACE_IF(Error <= ETol)(Error);

         //TRACE(Error);
         LinearAlgebra::Vector<complex> CoefficientsNew = X(LinearAlgebra::all, 0)[range(0,i)];
         Error = norm_2(direct_sum(Coefficients, complex(0.0)) - CoefficientsNew);
         //TRACE(Error)(CoefficientsNew);
         Coefficients = CoefficientsNew;
      }

      if (Beta[i] <= BetaTol)
      {
         DEBUG_TRACE("Beta hit tolerance")(Beta[i]);
         ++i;
         break;
      }
   }
   Iterations = i-1;

   // Calculate the matrix elements in the Krylov basis
#if 0
   LinearAlgebra::Matrix<double> M(Iterations, Iterations, 0.0);
   M(0,0) = Alpha[1];
   for (int k = 1; k < Iterations; ++k)
   {
      M(k,k) = Alpha[k+1];
      M(k,k-1) = M(k-1,k) = Beta[k];
   }
#endif

   // Exponentiate
   LinearAlgebra::Matrix<std::complex<double> > P = Theta * M(range(0,Iterations), range(0,Iterations));
   LinearAlgebra::Matrix<std::complex<double> > X = LinearAlgebra::Exponentiate(1.0, P);

   // expand the result vector
   VectorType Out = X(0,0) * v[1];
   for (unsigned k = 1; k < size1(X); ++k)
   {
      DEBUG_TRACE(X(k,0));
      Out += X(k,0) * v[k+1];
   }

   // Normalize
   Out *= Beta[0];
   DEBUG_TRACE(inner_prod(Out, Out));
   DEBUG_TRACE(inner_prod(x, x));

   ETol = Error;

   return Out;
}

// default value for ETol
template <typename VectorType, typename MultiplyFunctor>
inline
VectorType LanczosExponential(VectorType const& x,
                              MultiplyFunctor MatVecMultiply, int& Iterations,
                              std::complex<double> const& Theta)
{
   double ETol = 0.0;
   return LanczosExponential(x, MatVecMultiply, Iterations, Theta, ETol);
}

// A version of LanczosExponential where the first Krylov vector is already known
template <typename VectorType, typename MultiplyFunctor>
VectorType LanczosExponential(VectorType const& x, VectorType const& y,
                        MultiplyFunctor MatVecMultiply, int& Iterations,
                        std::complex<double> const& Theta)
{
   std::vector<VectorType>     v(Iterations+1);         // the Krylov vectors
   std::vector<double>         Alpha(Iterations+1);
   std::vector<double>         Beta(Iterations+1);

   VectorType r = x;
   Beta[0] = norm_frob(r);

   v[1] = (1.0 / Beta[0]) * r;
   r = y;
   Alpha[1] = real(inner_prod(v[1], r));
   r = r - Alpha[1] * v[1];
   Beta[1] = norm_frob(r);

   v[2] = (1.0 / Beta[1]) * r;
   r = MatVecMultiply(v[2]);
   r = r - inner_prod(v[1], r) * v[1];
   Alpha[2] = real(inner_prod(v[2], r));
   r = r - Alpha[2] * v[2];
   Beta[2] = norm_frob(r);

   double tol = 1E-14;

   int i;
   for (i = 3; i <= Iterations; ++i)
   {
      v[i] = (1.0 / Beta[i-1]) * r;
      r = MatVecMultiply(v[i]);
      r = r - Beta[i-1] * v[i-1];
      Alpha[i] = real(inner_prod(v[i], r));
      r = r - Alpha[i] * v[i];
      Beta[i] = norm_frob(r);

      if (Beta[i] <= tol)
      {
         DEBUG_TRACE("Beta hit tolerance")(Beta[i]);
         ++i;
         break;
      }
   }
   Iterations = i-1;

   // Calculate the matrix elements in the Krylov basis
   LinearAlgebra::Matrix<double> M(Iterations, Iterations, 0.0);
   M(0,0) = Alpha[1];
   for (int k = 1; k < Iterations; ++k)
   {
      M(k,k) = Alpha[k+1];
      M(k,k-1) = M(k-1,k) = Beta[k];
   }

   // Exponentiate
   LinearAlgebra::Matrix<std::complex<double> > P = Theta * M;
   LinearAlgebra::Matrix<std::complex<double> > X = LinearAlgebra::Exponentiate(1.0, P);

   // expand the result vector
   VectorType Out = X(0,0) * v[1];
   for (unsigned k = 1; k < size1(X); ++k)
   {
      DEBUG_TRACE(X(k,0));
      Out += X(k,0) * v[k+1];
   }

   // Normalize
   Out *= Beta[0];
   DEBUG_TRACE(inner_prod(Out, Out));
   DEBUG_TRACE(inner_prod(x, x));

   return Out;
}
