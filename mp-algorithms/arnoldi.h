// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/arnoldi.h
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
/*
   Arnoldi solver

   The preconditioning step is missing.

   A typical value for Tol is 1E-12, typical subspace size (iterations) is 20 - 30.

   If the iterations are stopped before convergence, the Tol on exit is set to
   the negative of the current residual norm.

   The SolverMode selects which eigenvalue we want, out of:
   LargestAlgebraicReal:  eigenvalue e with highest e.real (top of the spectrum)
   SmallestAlgebraicReal: eigenvalue e with lowest e.real  (bottom of the spectrum)
   LargestMagnitude:      eigenvalue e with largest |e|
   SmallestMagnitude:     eigenvalue e with smallest |e| (unstable unless the operator is positive)

   SmallestAlgebraicReal will find the same eigenvector that
   LargestAlgebraicReal would find if we used the negative of the operator

   TODO: we could forward this to ARPACK, if its available

*/

#include "common/proccontrol.h"
//#include "gram-schmidt.h"
#include <iostream>
#include <cmath>
#include "blas/functors.h"
#include "blas/matrix-eigen.h"

#if defined(SOLVER_TRACE_DETAILED)
#define TRACE_ARNOLDI(Msg) TRACE(Msg)
#else
#define TRACE_ARNOLDI(Msg) DUMMY_TRACE(Msg)
#endif

namespace LinearSolvers
{

double const DGKS_Threshold = 1.0 / std::sqrt(2.0); // 1.0; // between 0 and 1.

double const ArnoldiBetaTol = 1E-14;

enum SolverMode { LargestAlgebraicReal, LargestMagnitudeReal, LargestMagnitude, SmallestMagnitude };

template <typename VectorType, typename MultiplyFunctor>
std::complex<double> Arnoldi(VectorType& Guess, MultiplyFunctor MatVecMultiply, int& Iterations,
                             double& Tol, SolverMode Mode, bool Normalize = true, int Verbose = 0)
{
   typedef std::complex<double> complex;
   std::vector<VectorType> v;                                 // the subspace vectors
   std::vector<VectorType> Hv;                                // H * the subspace vectors
   blas::Matrix<complex>   SubH(Iterations, Iterations, 0.0); // matrix elements of H in the subspace
   complex Theta;         // eigenvalue
   v.reserve(Iterations);
   Hv.reserve(Iterations);

   VectorType w = copy(Guess);

   double Beta = norm_frob(w);
   double OrigBeta = Beta;
   w *= 1.0 / Beta;
   v.push_back(std::move(w));

   if (Verbose > 1)
   {
      std::cerr << "arnoldi: starting matrix-vector multiply\n";
      double Start = ProcControl::GetCPUTime();
      w = MatVecMultiply(v[0]);
      double CPU = ProcControl::GetCPUTime() - Start;
      std::cerr << "arnoldi: matrix-vector multiply took " << CPU << " seconds\n";
   }
   else
   {
      w = MatVecMultiply(v[0]);
   }
   SubH(0,0) = inner_prod(v[0], w);
   Hv.push_back(copy(w));
   w -= SubH(0,0) * v[0];

   Beta = norm_frob(w);
   TRACE_ARNOLDI(Beta);
   if (Beta < ArnoldiBetaTol)
   {
      if (Verbose > 0)
         std::cerr << "arnoldi: immediate return, invariant subspace found, Beta=" << Beta << '\n';
      TRACE_ARNOLDI("Immediate return - invariant subspace found")(Beta);
      Guess = std::move(v[0]);
      if (Normalize)
         Guess *= OrigBeta;
      Iterations = 1;
      Tol = Beta;
      return SubH(0,0);
   }

   for (int j = 1; j < Iterations; ++j)
   {
      SubH(j, j-1) = Beta;
      w *= 1.0 / Beta;
      v.push_back(std::move(w));

      // Matrix vector multiply
      if (Verbose > 1)
      {
         std::cerr << "arnoldi: starting matrix-vector multiply\n";
         double Start = ProcControl::GetCPUTime();
         w = MatVecMultiply(v[j]);
         double CPU = ProcControl::GetCPUTime() - Start;
         std::cerr << "arnoldi: matrix-vector multiply took " << CPU << " seconds\n";
      }
      else
      {
         w = MatVecMultiply(v[j]);
      }
      Hv.push_back(copy(w));
      // Subspace matrix elements
      double NormFrobSqH = 0;
      for (int i = 0; i <= j; ++i)
      {
         complex z = inner_prod(v[i], w);
	 //TRACE_ARNOLDI(z);
         SubH(i,j) = z;
         NormFrobSqH += norm_frob_sq(z);
         w -= z * v[i];
      }

      // DGKS correction
      double NormFrobSqF = norm_frob_sq(w);
      if (NormFrobSqF < DGKS_Threshold * DGKS_Threshold * NormFrobSqH)
      {
         TRACE_ARNOLDI("DGKS")(NormFrobSqF)(NormFrobSqH); //(SubH(range(0,j+1),range(0,j+1)));
         NormFrobSqH = 0;
         for (int i = 0; i <= j; ++i)
         {
            complex z = inner_prod(v[i], w);
	    //TRACE_ARNOLDI(z);
            SubH(i,j) += z;
            NormFrobSqH += norm_frob_sq(SubH(i,j));
            w -= z * v[i];
         }
         NormFrobSqF = norm_frob_sq(w);
	 TRACE_ARNOLDI("DGKS finished")(NormFrobSqF);
#if 0
         // attempt to detect breakdown of orthogonality - doesn't really work
         if (NormFrobSqF < DGKS_Threshold * DGKS_Threshold * NormFrobSqH)
         {
            // breakdown
            TRACE("Orthogonality Breakdown")(NormFrobSqF)(NormFrobSqH)(SubH(range(0,j+1),range(0,j+1)));
            for (int i = 0; i <= j; ++i)
            {
               TRACE(inner_prod(v[i], w));
            }
         }
#endif

      }

      blas::Matrix<complex> sH = SubH(blas::range(0,j+1), blas::range(0,j+1));
      //TRACE_ARNOLDI(sH);
      blas::Matrix<complex> Left(j,j), Right(j,j); // left and right eigenvectors
      blas::Vector<complex> Eigen(j);
      Diagonalize(std::move(sH), Eigen, Left, Right);

      //      TRACE(Eigen);
      int EigenIndex = 0;
      Theta = Eigen[EigenIndex];
      double ThetaMag = 0;
      switch (Mode)
      {
         case LargestMagnitudeReal : ThetaMag = norm_frob(Theta.real()); break;
         case LargestAlgebraicReal : ThetaMag = Theta.real(); break;
         case LargestMagnitude : ThetaMag = norm_frob(Theta); break;
         case SmallestMagnitude : ThetaMag = -norm_frob(Theta); break;
            //         case ClosestUnity : ThetaMag = -norm_frob(1.0 - Theta); break;
      }

      for (unsigned i = 1; i < Eigen.size(); ++i)
      {
         double NextMag = 0;
         switch (Mode)
         {
            case LargestMagnitudeReal : NextMag = norm_frob(Eigen[i].real()); break;
            case LargestAlgebraicReal : NextMag = Eigen[i].real(); break;
            case LargestMagnitude : NextMag = norm_frob(Eigen[i]); break;
            case SmallestMagnitude : NextMag = -norm_frob(Eigen[i]); break;
               //            case ClosestUnity : NextMag = -norm_frob(1.0 - Eigen[i]); break;
         }
         if (NextMag > ThetaMag)
         {
            EigenIndex = i;
            Theta = Eigen[i];
            ThetaMag = NextMag;
         }
      }

      TRACE_ARNOLDI(Eigen);

      // Calculate y = Ritz vector of the eigenvector
      VectorType y = Right(EigenIndex,0) * v[0];
      for (int i = 1; i <= j; ++i)
         y += Right(EigenIndex,i) * v[i];

      // Calculate the residual vector r = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int i = 0; i <= j; ++i)
         r += Right(EigenIndex,i) * Hv[i];

      double ResidNorm = norm_frob(r) / norm_frob(Theta);

      if (Verbose > 1)
         std::cerr << "arnoldi: iterations=" << (j+1)
                   << ", ResidNorm=" << ResidNorm
                   << ", evalue=" << Theta << '\n';

      TRACE_ARNOLDI(ResidNorm);

      if (ResidNorm < Tol)
      {
         if (Verbose > 0)
            std::cerr << "arnoldi: early return, residual norm below threshold, ResidNorm=" << ResidNorm
                      << ", iterations=" << (j+1) << '\n';
         TRACE_ARNOLDI("Early return - residual below threshold")(ResidNorm)(j);
         //      TRACE_ARNOLDI(Eigen);
         Guess = std::move(y);
         if (Normalize)
            Guess *= OrigBeta;
         Iterations = j+1;
         Tol = ResidNorm;
         return Theta;
      }
      //TRACE(norm_frob_sq(r));

      if (j == Iterations-1)  // finished?
      {
         if (Verbose > 0)
            std::cerr << "arnoldi: reached the maximum number of iterations, ResidNorm="
                      << ResidNorm << '\n';
         Guess = std::move(y);
         if (Normalize)
            Guess *= OrigBeta;
         //TRACE_ARNOLDI(Eigen);
         Tol = -ResidNorm;
         return Theta;
      }

      Beta = norm_frob(w);
      TRACE_ARNOLDI(Beta);
      if (Beta < ArnoldiBetaTol)
      {
         if (Verbose > 0)
            std::cerr << "arnoldi: early return, invariant subspace found, Beta=" << Beta
                      << ", iterations=" << (j+1) << '\n';
         TRACE_ARNOLDI("Early return - invariant subspace found")(Beta)(j);
         //      TRACE_ARNOLDI(Eigen);
         Guess = std::move(y);
         if (Normalize)
            Guess *= OrigBeta;
         Iterations = j+1;
         Tol = Beta;
         return Theta;
      }

   }

   return -1; // we never get here
}

// backwards-compatible version
template <typename VectorType, typename MultiplyFunctor>
inline
std::complex<double> Arnoldi(VectorType& Guess, MultiplyFunctor MatVecMultiply, int& Iterations)
{
   double Tol = 1E-14;
   return Arnoldi(Guess, MatVecMultiply, Iterations, Tol);
}

} // namespace LinearSolvers
