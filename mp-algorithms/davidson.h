// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/davidson.h
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
   Davidson solver.

   This can operate in 3 modes, which determines how we select the desired
   eigenvalue.
   DavidsonMode::Lowest is the 'standard' mode, which selects the arithmetic
   owest eigenvalue.
   DavidsonMode::Target selects the eigenvalue closest to the selected target
   eigenvalue.
   DavidsonMode::MaxOverlap selects the eigenvector that has the largest overlap
   with the initial guess vector.
*/

#include "linearalgebra/eigen.h"
#include "common/proccontrol.h"
#include <iostream>
#include <cmath>

// a hack for the preconditioner, should be moved somewhere else
#include "mps/state_component.h"

namespace LinearSolvers
{
//
using namespace LinearAlgebra;
typedef std::complex<double> complex;

enum class DavidsonMode { Lowest, Target, MaxOverlap };

template <typename VectorType>
bool GramSchmidtAppend(std::vector<VectorType>& Basis,
                       VectorType NewVec,
                       double Ortho = 2.0, double ParallelThreshold = 1E-32)
{
   double OriginalNorm = norm_frob_sq(NewVec);
   if (OriginalNorm < 1E-20)
      return false;

   int BasisSize = Basis.size();
   Basis.push_back(NewVec);

   double Norm2 = OriginalNorm;
   for (int n = 0; n < 2; ++n)
   {
      double MaxOverlap = 0.0;
      for (int i = 0; i < BasisSize; ++i)
      {
         complex Overlap = inner_prod(Basis[i], Basis.back());
         MaxOverlap = std::max(norm_frob_sq(Overlap), MaxOverlap);
         Basis.back() -= Overlap * Basis[i];
      }
      Norm2 = norm_frob_sq(Basis.back());

      if (Norm2 / OriginalNorm <= ParallelThreshold)  // parallel - cannot add the vector.
      {
         Basis.pop_back();
         return false;
      }
   }
   Basis.back() *= (1.0 / sqrt(Norm2));  // normalize
   return true;
}

struct InverseDiagonalPrecondition
{
   InverseDiagonalPrecondition(StateComponent const& Diag_, std::complex<double> Energy_)
      : Diag(Diag_), Energy(Energy_)
   {
      //      TRACE(Diag);
      for (unsigned i = 0; i < Diag.size(); ++i)
      {
	 for (StateComponent::operator_type::iterator I = iterate(Diag[i]); I;  ++I)
	 {
	    for (StateComponent::operator_type::inner_iterator J = iterate(I); J; ++J)
	    {
	       for (auto II = iterate(*J); II; ++II)
	       {
		  for (auto JJ = iterate(II); JJ; ++JJ)
		  {
                     // regularized inverse
                     *JJ = conj(*JJ - Energy) / (norm_frob_sq(*JJ-Energy) + 1E-8);
		  }
	       }
	    }
	 }
      }
      //TRACE(Diag)(Energy);
   }

   StateComponent operator()(StateComponent const& x) const
   {
      StateComponent Result(x);
      for (unsigned i = 0; i < x.size(); ++i)
      {
	 for (StateComponent::operator_type::iterator I = iterate(Result[i]); I;  ++I)
	 {
	    for (StateComponent::operator_type::inner_iterator J = iterate(I); J; ++J)
	    {
	       StateComponent::operator_type::const_inner_iterator d = iterate_at(Diag[i].data(), J.index1(), J.index2());
	       if (d)
	       {
		  *J = transform(*J, *d, LinearAlgebra::Multiplication<std::complex<double>, std::complex<double>>());
	       }
	    }
	 }
      }
      return Result;
   }

   StateComponent Diag;
   std::complex<double> Energy;
};


template <typename VectorType>
void
JacobiDavidsonIteration(VectorType& r, VectorType const& Diagonal, double Theta, VectorType const& Eigenvector)
{
   // orthogonalize against the Eigenvector
   r -= inner_prod(Eigenvector, r) * Eigenvector;
   // approximately solve (A-Theta)^{-1} x = r
   r = InverseDiagonalPrecondition(Diagonal, Theta)(r);
   // orthogonalize again
   r -= inner_prod(Eigenvector, r) * Eigenvector;
}

template <typename VectorType, typename MultiplyFunctor>
double Davidson(VectorType& Guess, VectorType const& Diagonal, MultiplyFunctor MatVecMultiply,
		DavidsonMode Mode, int& Iterations, int Verbose = 0, double TargetEnergy = 0.0,
                double MaxEnergy = 0.0)
{
   std::vector<VectorType> v;                                 // the subspace vectors
   std::vector<VectorType> Hv;                                // H * the subspace vectors
   Matrix<complex>         SubH(Iterations, Iterations, 0.0); // matrix elements of H in the subspace
   double Theta;         // eigenvalue
   v.reserve(Iterations);
   Hv.reserve(Iterations);

   MaxEnergy = TargetEnergy+1.0;

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
      //      TRACE(sH);
      Vector<double> Eigen = DiagonalizeHermitian(sH);
      //TRACE(sH)(Eigen);
      // The eigenvalues are ordered, so lowest energy is in Eigen[0]
      int n = 0;
      if (Mode == DavidsonMode::MaxOverlap)
      {
         // Find the vector with maximum overlap with the initial state (which is v[0]),
         // so this will be the vector n that has a maximum entry of sH(n,0)
         double MaxElement = norm_frob_sq(sH(n,0));
         double E = Eigen[n];
         for (int i = 1; i < j; ++i)
         {
            if ((E < TargetEnergy) || (Eigen[i] >= TargetEnergy && Eigen[i] <= MaxEnergy))
            {
               double ThisMaxElement = norm_frob_sq(sH(i,0));
               if (E < TargetEnergy || ThisMaxElement > MaxElement)
               {
                  MaxElement = ThisMaxElement;
                  n = i;
                  E = Eigen[n];
               }
            }
         }
         if (Verbose > 2)
            std::cerr << "Overlap " << MaxElement << '\n';
      }
      else if (Mode == DavidsonMode::Target)
      {
         double EDiff = norm_frob(Eigen[0] - TargetEnergy);
         for (int i = 1; i < j; ++i)
         {
            double ThisEDiff = norm_frob(Eigen[i] - TargetEnergy);
            if (ThisEDiff < EDiff)
            {
               EDiff = ThisEDiff;
               n = i;
            }
         }
      }
      // the eigenpair we want is Eigen[n], sH(n,all)
      Theta = Eigen[n];

      // Calculate y = Ritz vector of the eigenvector
      VectorType y = sH(n,0) * v[0];
      for (int i = 1; i < j; ++i)
         y += sH(n,i) * v[i];

      if (j == Iterations)  // finished?
      {
         Guess = y;
         Guess *= Beta; // normalize to the norm of the initial guess vector
         return Theta;
      }

      // Residual r = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int i = 0; i < j; ++i)
         r += sH(n,i) * Hv[i];

      // Insert preconditioning step here
      // Jacobi-Davidson would orthogonalize r against the Ritz eigenvector y here
      JacobiDavidsonIteration(r, Diagonal, Theta, y);

      // Orthogonalization step
      bool Added = GramSchmidtAppend(v, r, 1E-6);
      if (!Added)
      {
         if (Verbose > 1)
            std::cerr << "Happy breakdown: failed to add subspace vector\n";
         Guess = y;
         Guess *= Beta; // normalize to the norm of the initial guess vector
         Iterations = j;
         return Theta;
      }
   }

   return -1; // we never get here
}

} // namespace LinearSolvers
