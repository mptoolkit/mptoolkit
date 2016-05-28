// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/lanczos_old.h
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

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

using LinearAlgebra::range;

#define FAST_B

template <typename VectorType, typename MultiplyFunctor>
double Lanczos(VectorType& Guess, MultiplyFunctor MatVecMultiply, int MaxIterations)
{
   std::vector<VectorType>     f(MaxIterations+1);
   std::vector<double>         a(MaxIterations+1);
   std::vector<double>         b(MaxIterations+1);
   std::vector<double>         NormSq(MaxIterations+1);

   f[0] = Guess;
   double PsiNorm = norm_frob(Guess);
   std::cout.precision(14);

   // do n=0 iteration
   int n = 0;
   VectorType Hfn = MatVecMultiply(f[n]);
   NormSq[n] = norm_frob_sq(f[n]);
   a[n] = real(inner_prod(Hfn, f[n])) / NormSq[n];
   b[n] = 0;
   //   double Overlap = 1;
   double Energy = a[0];

   //   TRACE(Guess);
   //   TRACE(a[0])(NormSq[0]);

   if (MaxIterations == 0) return Energy;

   //   std::cout << "   n                    a                  b^2               Energy     Norm^2    Overlap\n";

   //      std::cout << std::setw(4) << n << " "
   //	<< std::setprecision(14)
   //	<< std::setw(20) << a[n] << " " << std::setw(20) << b[n]
   //	<< " " << std::setw(20) << Energy
   //	<< std::setprecision(4)
   //	<< " " << std::setw(10) << std::sqrt(NormSq[n]) << " " << std::setw(10)
   //	<< Overlap << std::endl;

   ++n; // n=1 iteration
   f[n] = Hfn - a[n-1] * f[n-1];
   NormSq[n] = norm_frob_sq(f[n]);

   Hfn = MatVecMultiply(f[n]);
   a[n] = real(inner_prod(Hfn, f[n])) / NormSq[n];
   b[n] = NormSq[n] / NormSq[n-1];
   //   Overlap = LinearAlgebra::norm_2(inner_prod(Guess, f[n])) / (PsiNorm * std::sqrt(NormSq[n]));

   LinearAlgebra::Matrix<std::complex<double> > M(MaxIterations+1, MaxIterations+1, 0.0);
   M(0,0) = a[0];
   M(1,1) = a[1];
   M(1,0) = M(0,1) = std::sqrt(b[1]);
   double EOld = Energy;
   //   TRACE(M(range(0,2), range(0,2)));
   Energy = LinearAlgebra::EigenvaluesHermitian(M(range(0,2), range(0,2)))[0];

   //     std::cout << std::setw(4) << n << " "
   //	<< std::setprecision(14)
   //	<< std::setw(20) << a[n] << " " << std::setw(20) << b[n]
   //	<< " " << std::setw(20) << Energy
   //	<< std::setprecision(4)
   //	<< " " << std::setw(10) << std::sqrt(NormSq[n]) << " " << std::setw(10)
   //	<< Overlap << std::endl;

   if (MaxIterations > 1)
   {
      while (n < MaxIterations && std::abs(Energy - EOld) > 1E-10)
      {
         ++n;
	 f[n] = Hfn - a[n-1] * f[n-1] - b[n-1] * f[n-2];
	 Hfn = MatVecMultiply(f[n]);
	 NormSq[n] = norm_frob_sq(f[n]);
	 a[n] = real(inner_prod(Hfn, f[n])) / NormSq[n];
	 b[n] = NormSq[n] / NormSq[n-1];
         //	 Overlap = LinearAlgebra::norm_2(inner_prod(Guess, f[n])) 
         //            / (PsiNorm * std::sqrt(NormSq[n]));

	 M(n,n) = a[n];
	 M(n-1,n) = M(n,n-1) = std::sqrt(b[n]);

	 EOld = Energy;
         //         TRACE(M(range(0,n+1), range(0,n+1)));
	 Energy = LinearAlgebra::EigenvaluesHermitian(M(range(0,n+1), range(0,n+1)))[0];

         // TRACE(Energy)(a[n])(inner_prod(Hfn, f[n]));
	 //	 double Next = ProcControl::GetCPUTime();
	 //	 Now = Next;
      }
   }

   TRACE(M);

   LinearAlgebra::Vector<double> EValues;
   EValues = DiagonalizeHermitian(M(range(0,n+1), range(0,n+1)));
   //   std::cout.precision(14);
   //   std::cout << "Lanczos eigenvalues are: " << EValues << std::endl;

   TRACE(EValues);

   // calculate the ground-state Lanczos vector
   Guess = (M(0,0) / std::sqrt(NormSq[0])) * f[0];
   for (int i = 1; i <= n && std::sqrt(NormSq[i]) > 1E-15; ++i)
   {
      //      TRACE(M(0,i))(NormSq[i]);
      Guess = Guess + (M(0,i) / std::sqrt(NormSq[i])) * f[i];
   }
   //   if (n > 1)
   //   {
   //      int i = n-1;
      //      TRACE(M(0,i))(f[i]);
   //      Guess = Guess + (M(0,i) / std::sqrt(NormSq[i])) * f[i];
   //   }

   // TRACE(norm_frob_sq(Guess));
   //   TRACE(Guess);

   return Energy;
}
