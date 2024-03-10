// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/davidson-extended.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
/*
   Davidson solver.  This is extended by allowing for an additional set of vectors
   that are outside the Hilbert space.  For these vectors, we must supply
   the matrix elements <i|j> and <i|H|j>.  For the subspace vectors |k>, we
   obtain <k|i> and <k|H|i>.  The Davidson algorithm then gives us the
   matrix elements <k1|H|k2>, with <k1|k2> = delta_{k1,k_2}.  A Cholesky
   decomposition then gives us an orthogonal matrix from which to determine
   required eigenpair.

   the OutVector' returns the output eigenvector as the coefficients of the extra
   vectors followed by the coefficient of the main component (returned in Guess').

   The preconditioning step is missing.
*/

#include "linearalgebra/eigen.h"
#include "linearalgebra/scalarmatrix.h"
#include "linearalgebra/matrix_utility.h"
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
   double OriginalNorm = norm_frob_sq(NewVec);
   if (OriginalNorm <= ParallelThreshold)
      return false;
   int BasisSize = Basis.size();
   Basis.push_back(NewVec);

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

Vector<double>
DiagonalizeNonOrtho(Matrix<complex>& Ident, Matrix<complex>& H)
{
   Matrix<complex> L = Ident;
   CholeskyFactorizeLower(L);
   zero_upper_triangular(L);
   Matrix<complex> Linv = L;
   InvertLowerTriangular(Linv);
   Matrix<complex> Ht = Linv * H * herm(Linv);
   Vector<double> Eigen = DiagonalizeHermitian(Ht);
   TRACE(Eigen);
   H = Ht * Linv;
   return Eigen;
}

template <typename VectorType, typename MultiplyFunctor>
double DavidsonExtended(VectorType& Guess, MultiplyFunctor MatVecMultiply,
                        std::vector<VectorType> const& ExtraVec,
                        std::vector<VectorType> const& ExtraHVec,
                        Matrix<complex> const& ExtraI,
                        Matrix<complex> const& ExtraH,
                        Vector<complex>& OutVector,
                        int& Iterations)
{
   int const ESize = ExtraVec.size();
   int const MaxSize = Iterations+ESize;  // Maximum size of the subspace
   std::vector<VectorType> v;                                 // the subspace vectors
   std::vector<VectorType> Hv;                                // H * the subspace vectors
   Matrix<complex>         SubH(MaxSize, MaxSize, 0.0); // matrix elements of H in the subspace
   Matrix<complex>         SubI(MaxSize, MaxSize, 0.0);  // <k|i>
   double Theta;         // eigenvalue

   v.reserve(Iterations);
   Hv.reserve(Iterations);

   // Initialize SubH and SubI
   SubI(range(0, ESize), range(0, ESize)) = ExtraI;
   SubI(range(ESize, MaxSize), range(ESize, MaxSize)) = ScalarMatrix<complex>(Iterations, Iterations, 1.0);
   SubH(range(0, ESize), range(0, ESize)) = ExtraH;

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
      // identity
      for (int i = 0; i < ESize; ++i)
      {
         complex z = inner_prod(ExtraVec[i], v[j-1]);
         SubI(i,j-1+ESize) = z;
         SubI(j-1+ESize,i) = conj(z);
         z = inner_prod(ExtraHVec[i], v[j-1]);
         SubH(i,j-1+ESize) = z;
         SubH(j-1+ESize,i) = conj(z);
      }
      for (int i = 0; i < j; ++i)
      {
         complex z = inner_prod(v[i], w);
         SubH(i+ESize,j-1+ESize) = z;
         SubH(j-1+ESize,i+ESize) = conj(z);
         w -= z * v[i];
      }

      // Subspace diagonalization
      Matrix<complex> sH = SubH(range(0,j+ESize), range(0,j+ESize));
      Matrix<complex> sI = SubI(range(0,j+ESize), range(0,j+ESize));
      Vector<double> Eigen = DiagonalizeNonOrtho(sI, sH);   // dominant eigenpair is Eigen[0],sH(0,all)
      Theta = Eigen[0];

      // Calculate y = Ritz vector of the eigenvector
      VectorType y = conj(sH(0,ESize)) * v[0];
      for (int i = 1; i < j; ++i)
      {
         y += conj(sH(0,ESize+i)) * v[i];
      }
      if (j == Iterations)  // finished?
      {
         Guess = y;
         double yNorm = norm_frob(y);  // this is less than 1 if we have sum vectors
         //Guess *= Beta; // normalize to the norm of the initial guess vector
         OutVector = Vector<complex>(ESize+1, 0.0);
         for (int i = 0; i < ESize; ++i)
         {
            OutVector[i] = conj(sH(0, i));
         }
         OutVector[ESize] = yNorm;
         return Theta;
      }

      // Residual r = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int i = 0; i < j; ++i)
         r += conj(sH(0,ESize+i)) * Hv[i];

      // Insert preconditioning step here

      // Orthogonalization step
      bool Added = GramSchmidtAppend(v, r, 1E-10);
      if (!Added)
      {
         WARNING("Failed to add subspace vector")(v.size());
         Guess = y;
         double yNorm = norm_frob(y);  // this is less than 1 if we have sum vectors
         //Guess *= Beta; // normalize to the norm of the initial guess vector
         Iterations = j;
         OutVector = Vector<complex>(ESize+1, 0.0);
         for (int i = 0; i < ESize; ++i)
         {
            OutVector[i] = conj(sH(0, i));
         }
         OutVector[ESize] = yNorm;
         return Theta;
      }
   }

   return -1; // we never get here
}

} // namespace LinearSolvers
