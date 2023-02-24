// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ef-matrix.h
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// Wrapper classes for the EA MPO solvers.
//
// These classes will calculate the E or F matrices for a triangular or string
// operator in the form of an InfiniteMPO.
//
// The ket and bra EA wavefunctions are stored in Psi(Tri)Upper/Lower
// respectively, with the infinite boundaries being stored in PsiUpper/Lower,
// and the window unit cells being stored in PsiTriUpper/Lower.
//
// The class calculates the (mixed) transfer matrix eigenvectors (TEVs) and E/F
// matrix elements as required and stores them in case they are needed later.
// Changing the operator or the windows (TODO) will automatically invalidate
// the effected matrix elements and eigenvalues.

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EF_MATRIX_H)
#define MPTOOLKIT_MP_ALGORITHMS_EF_MATRIX_H

#include "mp-algorithms/triangular_mpo_solver.h"
#include "mp-algorithms/transfer.h"

// Struct to hold the settings to initialize the EFMatrix class.
struct EFMatrixSettings
{
   int Degree = 0;
   double Tol = DefaultTol;
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   bool NeedFinalMatrix = true;
   bool EAOptimization = false;
   int Verbose = 0;
};

class EFMatrix
{
   public:
      EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings);

      // Set boundary wavefunctions: at the moment, we assume that the boundary
      // wavefunctions for the bra and ket are the same, as this allows us to
      // obtain the TEVs from the lambda matrix. (i >= 0)
      void SetPsi(int i, InfiniteWavefunctionLeft const& Psi);
      void SetPsi(int i, InfiniteWavefunctionRight const& Psi);

      // Set the window using a deque of single-site windows for each position
      // in the unit cell. (i > 0)
      void SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK);
      void SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK);

      // Update the operator: invalides the calculated E/F matrix elements and
      // (for string operators) the TEVs.
      void SetOp(InfiniteMPO Op_, int Degree_ = 0);

      // Returns the left and right (unit) TEVs, calculating them if they have
      // not been calculated yet. If the spectral radius is < 1, then
      // TLeft/Right will be null.
      MatrixOperator GetTLeft(int i, int j, bool F = false);
      MatrixOperator GetTRight(int i, int j, bool F = false);

      // Ident is the TEV in the same direction as the E/F matrix.
      MatrixOperator GetIdent(int i, int j, bool F = false)
         { return F ? this->GetTRight(i, j, true) : this->GetTLeft(i, j); }
      // Rho is the TEV in the opposite direction as the E/F matrix
      // (calculating the inner product of an E/F matrix element with Rho gives
      // the expectation value).
      MatrixOperator GetRho(int i, int j, bool F = false)
         { return F ? this->GetTLeft(i, j, true) : this->GetTRight(i, j); }

      // Get the E/F matrix where the ket and bra wavefunctions go up to the
      // boundaries i and j respectively.
      std::vector<KMatrixPolyType> GetE(int i, int j, int n);
      std::vector<KMatrixPolyType> GetE(int i, int j) { return this->GetE(i, j, -1); }
      std::vector<KMatrixPolyType> GetF(int i, int j, int n);
      std::vector<KMatrixPolyType> GetF(int i, int j) { return this->GetF(i, j, UnitCellSize); }

      std::vector<KMatrixPolyType> GetElement(int i, int j, bool F = false)
         { return F ? this->GetF(i, j) : this->GetE(i, j); }

   private:
      void CheckOperator();

      // Set the TEVs when PsiUpper and PsiLower are the same and the lambda
      // matrix for the left/right canonical form is known.
      void SetDiagTEVsLC(int i, RealDiagonalOperator Lambda);
      void SetDiagTEVsRC(int i, RealDiagonalOperator Lambda);

      // Solve for the the TEVs with an eigenvalue of magnitude 1: if the
      // spectral radius is less than one, set the TEVs to null.
      void CalculateTEVs(int i, int j);

      // Return the momentum factor corresponding to the element (i, j).
      std::complex<double> MomentumFactor(int i, int j) { return ExpIKUpper[i] * std::conj(ExpIKLower[j]); }

      // Calculate the unit cell for element (i, j).
      void CalculateE(int i, int j);
      void CalculateF(int i, int j);

      // Solve the corner element (i, j), using the off-diagonal component CTriK.
      void SolveE(int i, int j, std::vector<KMatrixPolyType> CTriK);
      void SolveF(int i, int j, std::vector<KMatrixPolyType> CTriK);

      std::map<int, LinearWavefunction> PsiUpper, PsiLower;
      std::map<int, std::map<int, StateComponent>> WindowUpper, WindowLower;
      std::map<int, std::complex<double>> ExpIKUpper, ExpIKLower;
      std::map<std::pair<int, int>, MatrixOperator> TLeft, TRight;
      std::map<std::pair<int, int>, bool> TCalculated;
      std::map<std::pair<int, int>, std::map<int, std::vector<KMatrixPolyType>>> EMatK, FMatK;

      int IMax, JMax;
      int UnitCellSize = 0;
      QuantumNumber QShift;

      InfiniteMPO Op;

      int Degree;
      double Tol;
      double UnityEpsilon;
      bool NeedFinalMatrix;
      bool EAOptimization;
      int Verbose;
};

#endif
