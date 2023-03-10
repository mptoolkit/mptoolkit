// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/triangular_mpo_solver_helpers.h
//
// Copyright (C) 2009-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TRIANGULAR_MPO_HELPERS_H)
#define MPTOOLKIT_MP_ALGORITHMS_TRIANGULAR_MPO_HELPERS_H

#include "wavefunction/operator_actions.h"
#include "mp-algorithms/arnoldi.h"
#include "mp-algorithms/gmres.h"
#include "common/environment.h"
#include "tensor/tensor_eigen.h"
#include "wavefunction/momentum_operations.h"

// Binomial coefficient
long Binomial(int n, int k);

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator const& Rhs, double Tol = 1E-14, int Verbose = 0)
{
   MatrixOperator Guess = MakeRandomMatrixOperator(Rhs.Basis1(), Rhs.Basis2(), Rhs.TransformsAs());
   LinearSolve(Guess, F, Rhs, Tol, Verbose);
   return Guess;
}

//#define USE_ITERATIVE_REFINEMENT
// iterative refinement doesn't appear to add anything useful here, but increasing the
// krylov length has a useful effect to avoid stagnation
template <typename Func>
void
LinearSolve(MatrixOperator& x, Func F, MatrixOperator const& Rhs, double Tol = 1E-14, int Verbose = 0)
{
   int const MaxIter = getenv_or_default("MP_GMRES_MAXITER", 10000);
   int m = 30;     // krylov subspace size
   int iter = 0;   // total number of iterations performed

   double normb = norm_frob(Rhs);

   int IterThisRound = m*50; // if it takes more than 50 rounds to converge, then we've proably stagnated
   double tol = Tol;
   int Ret = GmRes(x, F, normb, Rhs, m, IterThisRound, tol, LinearAlgebra::Identity<MatrixOperator>(), Verbose);
   iter += IterThisRound;

   while (Ret != 0 && iter < MaxIter)
   {
      // Attempt to avoid stagnation by increasing the number of iterations
      m=m+10;
      if (Verbose > 1)
      {
         std::cerr << "Refinement step, increasing m to " << m << '\n';
      }

      //      TRACE("Refinement step")(iter);
      // iterative refinement step
#if defined(USE_ITERATIVE_REFINEMENT)
      MatrixOperator R = Rhs- F(x);
      MatrixOperator xRefine = R;
#else
      MatrixOperator R = Rhs;
      MatrixOperator xRefine = x;
#endif
      IterThisRound = m*50;
      double tol = Tol;
      Ret = GmRes(xRefine, F, normb, R, m, IterThisRound, tol, LinearAlgebra::Identity<MatrixOperator>(), Verbose);

      iter += IterThisRound;
#if defined(USE_ITERATIVE_REFINEMENT)
      x += xRefine;
#else
      x = xRefine;
#endif

      if (Verbose > 1)
      {
         double Resid = norm_frob(F(x) - Rhs) / normb;
         std::cerr << "Residual after refinement step = " << Resid << '\n';
      }
   }

   if (Ret != 0)
   {
      // failed
      PANIC("Linear solver failed to converge after max_iter iterations")(MaxIter);
   }
}

// Calculate the (complex) eigenvalue that is closest to 1.0
// using Arnoldi.
// **DEPRECATED**
template <typename T>
std::complex<double>
FindClosestUnitEigenvalue(MatrixOperator& M, T Func, double tol, int Verbose)
{
   int Iterations = 20;
   double Tol = tol;
   std::complex<double> EtaL;
   EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, false, Verbose);
   while (Iterations == 20)
   {
      Tol = tol;
      EtaL = LinearSolvers::Arnoldi(M, Func, Iterations, Tol, LinearSolvers::LargestMagnitude, false, Verbose);
   }
   return EtaL;
}


// Given a matrix polynomial C, extract the components proportional to the
// UnitMatrixLeft into a c-number polynomial
ComplexPolyType
FindParallelParts(MatrixPolyType& C, MatrixOperator const& Identity, MatrixOperator const& Rho, int MaxDegree);

// Update the components of the E polynomial, given CParallel
ComplexPolyType
UpdateParallelParts(ComplexPolyType const& CParallel, std::complex<double> MomentumFactor = 1.0);

// Extract the parallel parts of C, and return the corresponding E
ComplexPolyType
DecomposeParallelParts(MatrixPolyType& C,
                       MatrixOperator const& Identity,
                       MatrixOperator const& Rho, double UnityEpsilon, int Degree);

// Extract the parallel parts of C, and return the corresponding E, with momentum.
// The Identity and Rho matrices can be any left/right eigenpair.
KComplexPolyType
DecomposeParallelPartsWithMomentum(KMatrixPolyType& C,
                                   std::complex<double> Factor,
                                   MatrixOperator const& Identity,
                                   MatrixOperator const& Rho, double UnityEpsilon, int Degree);

// Decompose the perpendicular parts at momentum K
// The Identity and Rho matrices can be any left/right eigenpair.
MatrixPolyType
DecomposePerpendicularPartsLeft(MatrixPolyType const& C, std::complex<double> K,
                               BasicFiniteMPO const& Diag,
                               MatrixOperator const& UnitMatrixLeft,
                               MatrixOperator const& UnitMatrixRight,
                               LinearWavefunction const& Psi1,
                               LinearWavefunction const& Psi2,
                               QuantumNumber const& QShift,
                               std::complex<double> Scale,
                               bool HasEigenvalue1,
                               double Tol,
                               int Verbose);

KMatrixPolyType
DecomposePerpendicularPartsLeft(KMatrixPolyType const& C,
                                BasicFiniteMPO const& Diag,
                                MatrixOperator const& UnitMatrixLeft,
                                MatrixOperator const& UnitMatrixRight,
                                LinearWavefunction const& Psi1,
                                LinearWavefunction const& Psi2,
                                QuantumNumber const& QShift,
                                std::complex<double> Scale,
                                bool HasEigenvalue1,
                                double Tol,
                                int Verbose);

// Decompose components perpendicular to the identity, in the right basis.
// Identity and Rho can be any right/left eigenpair
MatrixPolyType
DecomposePerpendicularPartsRight(MatrixPolyType const& C, std::complex<double> K,
                                 BasicFiniteMPO const& Diag,
                                 MatrixOperator const& Identity,
                                 MatrixOperator const& Rho,
                                 LinearWavefunction const& Psi1,
                                 LinearWavefunction const& Psi2,
                                 QuantumNumber const& QShift,
                                 std::complex<double> Scale,
                                 bool HasEigenvalue1,
                                 double Tol,
                                 int Verbose);

KMatrixPolyType
DecomposePerpendicularPartsRight(KMatrixPolyType const& C,
                                 BasicFiniteMPO const& Diag,
                                 MatrixOperator const& UnitMatrixLeft,
                                 MatrixOperator const& UnitMatrixRight,
                                 LinearWavefunction const& Psi1,
                                 LinearWavefunction const& Psi2,
                                 QuantumNumber const& QShift,
                                 std::complex<double> Scale,
                                 bool HasEigenvalue1,
                                 double Tol,
                                 int Verbose);

// Solve the components for the case where the diagonal operator is zero
KMatrixPolyType
SolveZeroDiagonal(KMatrixPolyType const& C);

// Solve components for the case where the diagonal operator is zero, at momentum K
MatrixPolyType
SolveZeroDiagonal(MatrixPolyType const& C, std::complex<double> K = 1.0);


struct SubProductLeftProject
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeftProject(LinearWavefunction const& Psi_, QuantumNumber const& QShift_,
                         MatrixOperator const& Proj_, MatrixOperator const& Ident_)
      : Psi(Psi_), QShift(QShift_), Proj(Proj_),
        Ident(Ident_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In; //delta_shift(In, QShift);
      for (LinearWavefunction::const_iterator I = Psi.begin(); I != Psi.end(); ++I)
       {
          Result = operator_prod(herm(*I), Result, *I);
       }
      Result = In - delta_shift(Result, QShift);
      //Result = 0.5 * (Result + adjoint(Result));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
};

struct SubProductRightProject
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRightProject(LinearWavefunction const& Psi_, QuantumNumber const& QShift_,
                          MatrixOperator const& Proj_, MatrixOperator const& Ident_)
      : Psi(Psi_), QShift(QShift_), Proj(Proj_), Ident(Ident_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In;
      LinearWavefunction::const_iterator I = Psi.end();
      while (I != Psi.begin())
      {
         --I;
         Result = operator_prod(*I, Result, herm(*I));
      }
      Result = delta_shift(Result, adjoint(QShift));
      Result = In - Result;
      //Result = 0.5 * (Result + adjoint(Result));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi;
   QuantumNumber QShift;
   MatrixOperator const& Proj;
   MatrixOperator const& Ident;
};

#endif
