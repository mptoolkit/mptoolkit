// -*- C++ -*- $Id$
/*
   Arnoldi solver

   The preconditioning step is missing.

   A typical value for Tol is 1E-12, typical subspace (iterations) is 20 - 30.

   If the iterations are performed without convergence, the Tol on exit is set to
   the negative of the current residual norm.

   The SolverMode selects which eigenvalue we want, out of:
   LargestAlgebraicReal:  eigenvalue e with highest e.real (top of the spectrum)
   SmallestAlgebraicReal: eigenvalue e with lowest e.real  (bottom of the spectrum)
   LargestMagnitude:      eigenvalue e with largest |e|

   SmallestAlgebraicReal will find the same eigenvector as LargestAlgebraicReal 
   with the negative matrix.

   TODO: we could forward this to ARPACK, if its available
   
*/

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
//#include "gram-schmidt.h"
#include <iostream>
#include <cmath>

#if defined(SOLVER_TRACE_DETAILED)
#define TRACE_ARNOLDI(Msg) TRACE(Msg)
#else
#define TRACE_ARNOLDI(Msg) DUMMY_TRACE(Msg)
#endif

namespace LinearSolvers
{

using namespace LinearAlgebra;

double const DGKS_Threshold = 1.0 / std::sqrt(2.0); // 1.0; // between 0 and 1.

double const ArnoldiBetaTol = 1E-14;

enum SolverMode { LargestAlgebraicReal, LargestMagnitudeReal, LargestMagnitude };

template <typename VectorType, typename MultiplyFunctor>
std::complex<double> Arnoldi(VectorType& Guess, MultiplyFunctor MatVecMultiply, int& Iterations,
			     double& Tol, SolverMode Mode, bool Normalize = true, bool Verbose = false)
{
   typedef std::complex<double> complex;
   std::vector<VectorType> v;                                 // the subspace vectors
   std::vector<VectorType> Hv;                                // H * the subspace vectors
   Matrix<complex>         SubH(Iterations, Iterations, 0.0); // matrix elements of H in the subspace
   complex Theta;         // eigenvalue
   v.reserve(Iterations);
   Hv.reserve(Iterations);
   
   VectorType w = Guess;

   double Beta = norm_frob(w);
   double OrigBeta = Beta;
   w *= 1.0 / Beta;
   v.push_back(w);

   w = MatVecMultiply(v[0]);
   SubH(0,0) = inner_prod(v[0], w);
   Hv.push_back(w);
   w -= SubH(0,0) * v[0];

   Beta = norm_frob(w);
   TRACE_ARNOLDI(Beta);
   if (Beta < ArnoldiBetaTol)
   {
      if (Verbose)
         std::cerr << "arnoldi: immediate return, invariant subspace found, Beta=" << Beta << '\n';
      TRACE_ARNOLDI("Immediate return - invariant subspace found")(Beta);
      Guess = v[0];
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
      v.push_back(w);

      // Matrix vector multiply
      w = MatVecMultiply(v[j]);
      Hv.push_back(w);
      // Subspace matrix elements
      double NormFrobSqH = 0;
      for (int i = 0; i <= j; ++i)
      {
	 complex z = inner_prod(v[i], w);
	 SubH(i,j) = z;
	 NormFrobSqH += LinearAlgebra::norm_frob_sq(z);
	 w -= z * v[i];
      }

      // DGKS correction
      double NormFrobSqF = norm_frob_sq(w);
      if (NormFrobSqF < DGKS_Threshold * DGKS_Threshold * NormFrobSqH)
      {
         TRACE_ARNOLDI("DGKS")(NormFrobSqF)(NormFrobSqH)(SubH(range(0,j+1),range(0,j+1)));
         NormFrobSqH = 0;
	 for (int i = 0; i <= j; ++i)
	 {
	    complex z = inner_prod(v[i], w);
	    SubH(i,j) += z;
            NormFrobSqH += LinearAlgebra::norm_frob_sq(SubH(i,j));
	    w -= z * v[i];
	 }
         NormFrobSqF = norm_frob_sq(w);

#if 0
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

      Matrix<complex> sH = SubH(range(0,j+1), range(0,j+1));
      //TRACE_ARNOLDI(sH);
      Matrix<complex> Left, Right; // left and right eigenvectors
      Vector<complex> Eigen = Diagonalize(sH, Left, Right);

      int EigenIndex = 0;
      Theta = Eigen[EigenIndex];
      double ThetaMag = 0;
      switch (Mode)
      {
         case LargestMagnitudeReal : ThetaMag = norm_frob(Theta.real()); break;
         case LargestAlgebraicReal : ThetaMag = Theta.real(); break;
         case LargestMagnitude : ThetaMag = norm_frob(Theta); break;
	    //         case ClosestUnity : ThetaMag = -norm_frob(1.0 - Theta); break;
      }

      for (unsigned i = 1; i < size(Eigen); ++i)
      {
	 double NextMag = 0;
	 switch (Mode)
	 {
            case LargestMagnitudeReal : NextMag = norm_frob(Eigen[i].real()); break;
            case LargestAlgebraicReal : NextMag = Eigen[i].real(); break;
   	    case LargestMagnitude : NextMag = norm_frob(Eigen[i]); break;
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

#if 1
      // Residual r = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int i = 0; i <= j; ++i)
	 r += Right(EigenIndex,i) * Hv[i];

      double ResidNorm = norm_frob(r);

      TRACE_ARNOLDI(ResidNorm);

      if (ResidNorm < Tol)
      {
         if (Verbose)
            std::cerr << "arnoldi: early return, residual norm below threshold, ResidNorm=" << ResidNorm 
                      << ", iterations=" << (j+1) << '\n';
	 TRACE_ARNOLDI("Early return - residual below threshold")(ResidNorm)(j);
	 //	 TRACE_ARNOLDI(Eigen);
	 Guess = y;
	 if (Normalize)
	    Guess *= OrigBeta;
	 Iterations = j+1;
         Tol = ResidNorm;
	 return Theta;
      }
      //TRACE(norm_frob_sq(r));
#endif

      if (j == Iterations-1)  // finished?
      {
	 Guess = y;
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
         if (Verbose)
            std::cerr << "arnoldi: early return, invariant subspace found, Beta=" << Beta 
                      << ", iterations=" << (j+1) << '\n';
	 TRACE_ARNOLDI("Early return - invariant subspace found")(Beta)(j);
	 //	 TRACE_ARNOLDI(Eigen);
	 Guess = y;
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
