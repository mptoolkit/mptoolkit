// -*- C++ -*- $Id$

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
	       double& Tol, int MinIter = 2, bool Verbose = false)
{
   std::vector<VectorType>     v;          // the Krylov vectors
   std::vector<VectorType>     Hv;         // the Krylov vectors

   LinearAlgebra::Matrix<double> SubH(Iterations+1, Iterations+1, 0.0);

   VectorType w = Guess;

   double Beta = norm_frob(w);
   CHECK(!isnan(Beta));
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
      SubH(i,i) = real(inner_prod(v[i], w));
      w -= SubH(i,i) * v[i];
      Beta = norm_frob(w);


      if (Beta < LanczosBetaTol)
      {
         // Early return, we can't improve over the previous energy and eigenvector
	 if (Verbose)
	    std::cerr << "lanczos: early return, invariant subspace found, Beta="
		      << Beta << ", iterations=" << (i+1) << '\n';
	 Iterations = i+1;
         LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i),
                                                LinearAlgebra::range(0,i));
         LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
         double Theta = EValues[0];    // smallest eigenvalue
         VectorType y = M(0,0) * v[0];
         for (int j = 1; j <= i; ++j)
            y += M(0,j) * v[j];
    	 Tol = Beta;
	 Guess = y;
	 return Theta;
      }

      // solution of the tridiagonal subproblem
      LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i+1),
					     LinearAlgebra::range(0,i+1));
      if (isnan(M(0,0)))
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
      VectorType y = M(0,0) * v[0];
      for (int j = 1; j <= i; ++j)
	 y += M(0,j) * v[j];

      // residual = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int j = 0; j < i; ++j)
	 r += M(0,j) * Hv[j];

      double ResidNorm = norm_frob(r);

      if (ResidNorm < fabs(Tol * Theta) && i > MinIter)
      {
	 if (Verbose)
	    std::cerr << "lanczos: early return, residual norm below threshold, ResidNorm="
		      << ResidNorm << ", iterations=" << (i+1) << '\n';
	 Iterations = i+1;
	 Tol = ResidNorm/fabs(Theta);
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

// backwards-compatible version
template <typename VectorType, typename MultiplyFunctor>
double Lanczos(VectorType& Guess, MultiplyFunctor MatVecMultiply, int Iterations)
{
   int Iter = Iterations;
   double Tol = 1E-10;
   return Lanczos(Guess, MatVecMultiply, Iter, Tol);
}

#endif
