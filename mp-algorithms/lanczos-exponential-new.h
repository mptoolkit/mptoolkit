// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/lanczos-exponential-new.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

template <typename VectorType, typename MultiplyFunctor>
VectorType LanczosExponential(VectorType const& x,
                              MultiplyFunctor MatVecMultiply, int& Iterations,
                              std::complex<double> const& Theta, double& ErrTol)
{
   typedef std::complex<double> complex;

   std::vector<VectorType> v(Iterations+1);
   std::vector<double> Alpha(Iterations+1);
   std::vector<double> Beta(Iterations);
   LinearAlgebra::Matrix<double> T(Iterations+1, Iterations+1, 0.0);
   LinearAlgebra::Vector<complex> C(1, 1.0);
   double Err = 1.0;
   double const BetaTol = 1e-14;

   v[0] = x * (1.0/norm_frob(x));

   VectorType w = MatVecMultiply(v[0]);
   Alpha[0] = real(inner_prod(w, v[0]));
   w = w - Alpha[0] * v[0];

   T(0, 0) = Alpha[0];

   int i;
   for (i = 1; i <= Iterations && Err > ErrTol; ++i)
   {
      Beta[i-1] = norm_frob(w);

      if (Beta[i-1] < BetaTol)
         break;

      v[i] = w * (1.0/Beta[i-1]);

      w = MatVecMultiply(v[i]);
      Alpha[i] = real(inner_prod(w, v[i]));
      w = w - Alpha[i] * v[i] - Beta[i-1] * v[i-1];

      T(i, i) = Alpha[i];
      T(i-1, i) = T(i, i-1) = Beta[i-1];

      LinearAlgebra::Matrix<complex> P = Theta * T(LinearAlgebra::range(0, i+1),
                                                   LinearAlgebra::range(0, i+1));
      LinearAlgebra::Matrix<complex> X = LinearAlgebra::Exponentiate(1.0, P);
      LinearAlgebra::Vector<complex> CNew = X(LinearAlgebra::all, 0);

      Err = norm_2(direct_sum(C, complex(0.0)) - CNew);

      C = CNew;
   }

   Iterations = i-1;
   ErrTol = Err;

   // Calculate the time-evolved vector.
   VectorType Out = v[0] * C[0];
   for (int j = 1; j <= Iterations; ++j)
   {
      Out += v[j] * C[j];
   }

   // Normalize Out.
   Out *= 1.0/norm_frob(Out);

   return Out;
}
