// -*- C++ -*- $Id$

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;

template <typename VectorType, typename MultiplyFunctor>
void LanczosExponential(VectorType const& x,
                        MultiplyFunctor MatVecMultiply, int Iterations,
                        sd::vector<std::complex<double> > const& Theta, 
                        std::vector<VectorType>& Out)
{
   std::vector<VectorType>     v(Iterations+1);         // the Krylov vectors
   std::vector<double>         Alpha(Iterations+1);
   std::vector<double>         Beta(Iterations+1);

   VectorType r = x;
   Beta[0] = norm_frob(r);

   double tol = 1E-14;

   int i;
   for (i = 1; i <= Iterations; ++i)
   {
      v[i] = (1.0 / Beta[i-1]) * r;
      r = MatVecMultiply(v[i]);
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

   // calculate the exponential
   Out.resize(Theta.size());
   for (int j = 0; j < Theta.size(); ++j)
   {
      std::complex<double> Th =  exp(Theta[j]) * Beta[0];
      Out[j] = (M(0,0) * Th) * v[1];
      for (int k = 1; k < Iterations; ++k)
      {
         Out[j] += (M(0,k) * Th) * v[k+1];
      }
   }
}
