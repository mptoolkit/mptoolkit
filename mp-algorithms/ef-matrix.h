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
//
// A note on the indices: the window of index i is between the boundaries i-1
// and i. The first boundary index is 0, and so the first window index is 1.
// For E matrices, the index signifies the position in the wavefunction from
// left to right, but for F matrices, it signifies the position from right to
// left.

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EF_MATRIX_H)
#define MPTOOLKIT_MP_ALGORITHMS_EF_MATRIX_H

#include "mp-algorithms/transfer.h"
#include "mp-algorithms/triangular_mpo_solver.h"

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

// Generic class to handle the code common to the E and F matrices.
class EFMatrix
{
   public:
      EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings);

      // Set boundary wavefunctions: at the moment, we assume that the boundary
      // wavefunctions for the bra and ket are the same, as this allows us to
      // obtain the TEVs from the lambda matrix. (i >= 0)
      void SetPsi(int i, InfiniteWavefunctionLeft const& Psi);
      void SetPsi(int i, InfiniteWavefunctionRight const& Psi);

      // Set the window unit cell between boundaries i-1 and i with the
      // (individual) momentum factor ExpIK. (i > 0)
      void SetPsiTriUpper(int i, LinearWavefunction const& PsiTri, std::complex<double> ExpIK);
      void SetPsiTriLower(int j, LinearWavefunction const& PsiTri, std::complex<double> ExpIK);

      // Set the window using a deque of single-site windows for each position
      // in the unit cell: requires the boundaries at i-1 and i to be specified.
      virtual void SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK) = 0;
      virtual void SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK) = 0;

      // Update the operator: invalides the calculated E/F matrix elements and
      // (for string operators) the TEVs.
      void SetOp(InfiniteMPO Op_, int Degree_ = 0);

      // Returns the left and right (unit) TEVs, calculating them if they have
      // not been calculated yet. If the spectral radius is < 1, then
      // TLeft/Right will be null.
      MatrixOperator GetTLeft(int i, int j);
      MatrixOperator GetTRight(int i, int j);

      // Ident is the TEV in the same direction as the E/F matrix.
      virtual MatrixOperator GetIdent(int i, int j) = 0;
      // Rho is the TEV in the opposite direction as the E/F matrix
      // (calculating the inner product of an E/F matrix element with Rho gives
      // the expectation value).
      virtual MatrixOperator GetRho(int i, int j) = 0;

      // Get the E/F matrix where the ket and bra wavefunctions go up to the
      // boundaries i and j respectively.
      std::vector<KMatrixPolyType> GetElement(int i, int j);

   protected:
      virtual void CheckOperator() = 0;

      virtual void SetDiagTEVsLC(int i, RealDiagonalOperator Lambda) = 0;
      virtual void SetDiagTEVsRC(int i, RealDiagonalOperator Lambda) = 0;

      virtual void CalculateTEVs(int i, int j) = 0;

      virtual void CalculateElement(int i, int j) = 0;

      std::map<int, LinearWavefunction> PsiUpper, PsiLower;
      std::map<int, LinearWavefunction> PsiTriUpper, PsiTriLower;
      std::map<int, std::complex<double>> ExpIKUpper, ExpIKLower;
      std::map<std::pair<int, int>, MatrixOperator> TLeft, TRight;
      std::map<std::pair<int, int>, bool> TCalculated;
      std::map<std::pair<int, int>, std::vector<KMatrixPolyType>> EFMatK;

      int IMax, JMax;
      QuantumNumber QShift;

      InfiniteMPO Op;

      int Degree;
      double Tol;
      double UnityEpsilon;
      bool NeedFinalMatrix;
      bool EAOptimization;
      int Verbose;
};

// The specific classes for E matrices and F matrices.
class EMatrix : public EFMatrix
{
   public:
      EMatrix(InfiniteMPO Op_, EFMatrixSettings Settings)
         : EFMatrix(Op_, Settings) { this->CheckOperator(); }

      void SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK);
      void SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK);

      MatrixOperator GetIdent(int i, int j) { return this->GetTLeft(i, j); }
      MatrixOperator GetRho(int i, int j) { return this->GetTRight(i, j); }

   protected:
      void CheckOperator();

      void SetDiagTEVsLC(int i, RealDiagonalOperator Lambda);
      void SetDiagTEVsRC(int i, RealDiagonalOperator Lambda);

      void CalculateTEVs(int i, int j);

      void CalculateElement(int i, int j);
};

class FMatrix : public EFMatrix
{
   public:
      FMatrix(InfiniteMPO Op_, EFMatrixSettings Settings)
         : EFMatrix(Op_, Settings) { this->CheckOperator(); }

      void SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK);
      void SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK);

      MatrixOperator GetIdent(int i, int j) { return this->GetTRight(i, j); }
      MatrixOperator GetRho(int i, int j) { return this->GetTLeft(i, j); }

   protected:
      void CheckOperator();

      void SetDiagTEVsLC(int i, RealDiagonalOperator Lambda);
      void SetDiagTEVsRC(int i, RealDiagonalOperator Lambda);

      void CalculateTEVs(int i, int j);

      void CalculateElement(int i, int j);
};

#endif
