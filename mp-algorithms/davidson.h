// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/davidson.h
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
/*
   Davidson solver

   The preconditioning step is missing.
*/

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

namespace LinearSolvers
{

using namespace LinearAlgebra;
typedef std::complex<double> complex;

template <typename VectorType>
bool GramSchmidtAppend(std::vector<VectorType>& Basis, 
		       VectorType NewVec, 
		       double Ortho = 2.0, double ParallelThreshold = 0.0)
{
   int BasisSize = Basis.size();
   Basis.push_back(NewVec);
   
   double OriginalNorm = norm_frob_sq(NewVec);
   double Norm2 = OriginalNorm;
   bool Converged = false;
   while (!Converged)
   {
      double MaxOverlap = 0.0;
      for (int i = 0; i < BasisSize; ++i)
      {
	 complex Overlap = inner_prod(Basis[i], Basis.back());
	 MaxOverlap = std::max(norm_frob_sq(Overlap), MaxOverlap);
         Basis.back() -= Overlap * Basis[i];
      }
      double NewNorm2 = norm_frob_sq(Basis.back());

      if (NewNorm2 / OriginalNorm <= ParallelThreshold)  // parallel - cannot add the vector.
      {
	 Basis.pop_back();
	 return false;
      }

      Converged = (MaxOverlap <= Ortho * sqrt(Norm2));
      Norm2 = NewNorm2;
   }
   Basis.back() *= (1.0 / sqrt(Norm2));  // normalize
   return true;
}

template <typename VectorType, typename MultiplyFunctor>
double Davidson(VectorType& Guess, VectorType const& Diagonal, MultiplyFunctor MatVecMultiply, int& Iterations)
{
   std::vector<VectorType> v;                                 // the subspace vectors
   std::vector<VectorType> Hv;                                // H * the subspace vectors
   Matrix<complex>         SubH(Iterations, Iterations, 0.0); // matrix elements of H in the subspace
   double Theta;         // eigenvalue
   v.reserve(Iterations);
   Hv.reserve(Iterations);
   
   VectorType w = Guess;

   double Beta = norm_frob(w);
   w *= 1.0 / Beta;
   v.push_back(w);

   for (int j = 1; j <= Iterations; ++j)
   {
      // Matrix vector multiply
      w = MatVecMultiply(v[j-1]);
      Hv.push_back(w);
      // Subspace matrix elements
      for (int i = 0; i < j; ++i)
      {
	 complex z = inner_prod(v[i], w);
	 SubH(i,j-1) = z;
	 SubH(j-1,i) = conj(z);
	 w -= z * v[i];
      }

      Matrix<complex> sH = SubH(range(0,j), range(0,j));
      Vector<double> Eigen = DiagonalizeHermitian(sH);   // dominant eigenpair is Eigen[0],sH(0,all)
      Theta = Eigen[0];

      // Calculate y = Ritz vector of the eigenvector
      VectorType y = sH(0,0) * v[0];
      for (int i = 1; i < j; ++i)
	 y += sH(0,i) * v[i];

      if (j == Iterations)  // finished?
      {
	 Guess = y;
	 Guess *= Beta; // normalize to the norm of the initial guess vector
	 return Theta;
      }
      
      // Residual r = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int i = 0; i < j; ++i)
	 r += sH(0,i) * Hv[i];

      // Insert preconditioning step here
      Precondition(r, Diagonal, Theta);

      // Orthogonalization step
      bool Added = GramSchmidtAppend(v, r, 1E-6);
      if (!Added)
      {
	 WARNING("Failed to add subspace vector")(v.size());
	 Guess = y;
	 Guess *= Beta; // normalize to the norm of the initial guess vector
	 Iterations = j;
	 return Theta;
      }
   }

   return -1; // we never get here
}

} // namespace LinearSolvers
