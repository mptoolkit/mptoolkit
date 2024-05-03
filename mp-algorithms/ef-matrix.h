// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
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
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

// Wrapper classes for the EA MPO solvers.
//
// These classes will calculate the E or F matrices for a triangular or string
// operator in the form of an InfiniteMPO.
//
// The ket and bra boundary wavefunctions are stored in PsiUpper/Lower
// respectively, and the momentum windows are stored in WindowUpper/Lower.
//
// The E/F matrix indices are represented using the nonnegative integers and
// "infinity" (which is implemented in the ExtendedInt class). The left and
// right boundaries have indices zero and infinity respectively, while the
// elements with nonzero, finite values are those with the window partially
// incorporated.
//
// Since the boundaries only exist for the index = zero or infiinty, we
// represent the index for objects which store the boundaries or boundary
// transfer matrices as bools, where false corresponds to zero and true
// corresponds to infinity (bool will implicity convert to ExtendedInt in this
// way).
//
// The indices for the windows are a bit tricky, since each window position
// corresponds to a FSM transition between two indices for the E/F matrix.
// The index for the first site in the window is 1, since this will avoid
// ambiguity once multiwindow excitations are implemented. The final index of
// the window is stored as I/JMax (and is not infinity), and the finite E/F
// matrix element indices only go up to I/JMax-1.
//
// The class calculates the (mixed) transfer matrix eigenvectors (TEVs) and E/F
// matrix elements as required and stores them in case they are needed later.
// Changing the operator or the windows will automatically invalidate the
// affected matrix elements and eigenvalues.

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EF_MATRIX_H)
#define MPTOOLKIT_MP_ALGORITHMS_EF_MATRIX_H

#include "mp-algorithms/triangular_mpo_solver.h"
#include "mp-algorithms/transfer.h"

// Struct to hold the settings to initialize the EFMatrix class.
struct EFMatrixSettings
{
   // Number of windows in the upper and lower states.
   int NUpper = 1;
   int NLower = 1;

   int Degree = 0;
   double Tol = DefaultTol;
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   bool NeedFinalMatrix = true;
   bool EAOptimization = false;
   bool SubtractEnergy = false;
   // Print a warning if we cannot automatically fix the phase of the TEVs.
   bool PhaseWarnings = true;
   int Verbose = 0;
};

// A class representing either a value of "zero" or "infinity" (essentially a
// wrapper for bool). This is used for represeting the indices for corner
// elements for EFMatrix.
class ZeroInf
{
   public:
      ZeroInf(bool Inf_)
         : Inf(Inf_) {}

      bool is_inf() const { return Inf; }

   private:
      bool Inf;
};

// A constant representing the value infinity (zero is just represented by "0",
// which can be implcitly convert to a ZeroInf object).
const ZeroInf Infinity(true);

// A class to represent an extended integer, which can be a finite integer or infinity.
class ExtendedInt
{
   public:
      ExtendedInt() {}

      ExtendedInt(int Value_)
         : Value(Value_), Inf(false) {}

      ExtendedInt(ZeroInf ZI)
         : Value(0), Inf(ZI.is_inf()) {}

      bool is_inf() const { return Inf; }

      // This should only be called if the number is finite.
      int value() const
      {
         if (Inf)
            throw std::runtime_error("ExtendedInt: called value for an infinite number.");
         return Value;
      }

   private:
      int Value;
      bool Inf;
};

class EFMatrix
{
   public:
      typedef std::vector<ExtendedInt> EAIndex;
      typedef std::vector<ZeroInf> CornerIndex;

      EFMatrix() {}

      EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings);

      // Set boundary wavefunctions: at the moment, we assume that the boundary
      // wavefunctions for the bra and ket are the same, as this allows us to
      // obtain the TEVs from the lambda matrix.
      void SetPsi(CornerIndex i, InfiniteWavefunctionLeft const& Psi, std::complex<double> ExpIK = 1.0);
      void SetPsi(CornerIndex i, InfiniteWavefunctionRight const& Psi, std::complex<double> ExpIK = 1.0);

      // Set the momentum factor.
      void SetExpIKUpper(CornerIndex i, std::complex<double> ExpIK);
      void SetExpIKLower(CornerIndex i, std::complex<double> ExpIK);

      // Set the windows using a vector of LinearWavefunctions.
      void SetWindowUpper(EAIndex i, std::vector<LinearWavefunction> const& WindowVec);
      void SetWindowLower(EAIndex j, std::vector<LinearWavefunction> const& WindowVec);

      // By default, set the diagonal window.
      void SetWindowUpper(std::vector<LinearWavefunction> const& WindowVec)
         { this->SetWindowUpper(EAIndex(NUpper, 1), WindowVec); }
      void SetWindowLower(std::vector<LinearWavefunction> const& WindowVec)
         { this->SetWindowLower(EAIndex(NLower, 1), WindowVec); }

      // Set a diagonal window element using a deque of single-site windows for
      // each position in the unit cell (in window order).
      void SetWUpper(int i, std::deque<StateComponent> const& BDeque);
      void SetWLower(int j, std::deque<StateComponent> const& BDeque);

      // Update the operator: invalidates the calculated E/F matrix elements and
      // (for string operators) the TEVs.
      void SetOp(InfiniteMPO Op_, int Degree_ = 0);

      // Returns the left and right (unit) TEVs, calculating them if they have
      // not been calculated yet. If the spectral radius is < 1, then
      // TLeft/Right will be null.
      MatrixOperator GetTLeft(CornerIndex i, CornerIndex j, bool F = false);
      MatrixOperator GetTRight(CornerIndex i, CornerIndex j, bool F = false);

      // Calculate all of the uncalculated transfer matrix eigenvalues: use
      // this if the EFMatrix class needs to be passed to ARPACK, since we need
      // ARPACK to calculate the TEVs and ARPACK isn't reentrant.
      void CalculateAllTEVs();

      // Ident is the TEV in the same direction as the E/F matrix.
      MatrixOperator GetIdent(CornerIndex i, CornerIndex j, bool F = false)
         { return F ? this->GetTRight(i, j, true) : this->GetTLeft(i, j); }
      // Rho is the TEV in the opposite direction as the E/F matrix
      // (calculating the inner product of an E/F matrix element with Rho gives
      // the expectation value).
      MatrixOperator GetRho(CornerIndex i, CornerIndex j, bool F = false)
         { return F ? this->GetTLeft(i, j, true) : this->GetTRight(i, j); }

      // Get the E/F matrix where the ket and bra wavefunctions go up to the
      // boundaries i and j respectively.
      std::vector<KMatrixPolyType> GetE(EAIndex i, EAIndex j, int n);
      std::vector<KMatrixPolyType> GetE(EAIndex i, EAIndex j) { return this->GetE(i, j, -1); }
      std::vector<KMatrixPolyType> GetF(EAIndex i, EAIndex j, int n);
      std::vector<KMatrixPolyType> GetF(EAIndex i, EAIndex j) { return this->GetF(i, j, UnitCellSize); }

      std::vector<KMatrixPolyType> GetElement(EAIndex i, EAIndex j, bool F = false)
         { return F ? this->GetF(i, j) : this->GetE(i, j); }

      // Get the constant, zero-momentum part of this element as a StateComponent.
      StateComponent GetESC(EAIndex i, EAIndex j, int n);
      StateComponent GetFSC(EAIndex i, EAIndex j, int n);

      // Calculate action of the effective Hamiltonian on the window w.r.t.
      // site i for each position in the unit cell (in window order).
      std::deque<StateComponent> GetHEff(int i = 0);

   private:
      void CheckOperator();

      // Set the TEVs when PsiUpper and PsiLower are the same and the lambda
      // matrix for the left/right canonical form is known.
      void SetDiagTEVsLC(CornerIndex i, RealDiagonalOperator Lambda);
      void SetDiagTEVsRC(CornerIndex i, RealDiagonalOperator Lambda);

      // Solve for the TEVs with an eigenvalue of magnitude 1: if the
      // spectral radius is less than one, set the TEVs to null.
      void CalculateTEVs(CornerIndex i, CornerIndex j);

      // Return the momentum factor corresponding to the element (i, j).
      std::complex<double> MomentumFactor(EAIndex i, EAIndex j);
      std::complex<double> MomentumFactor(CornerIndex i, CornerIndex j);

      // Return the StateComponent for the window at position n in the unit
      // cell for the index i/j.
      StateComponent GetWUpper(EAIndex i, int n);
      StateComponent GetWLower(EAIndex j, int n);
      StateComponent GetW(EAIndex i, int n, bool Upper)
         { return Upper ? this->GetWUpper(i, n) : this->GetWLower(i, n); }

      // Get a map of all possible indices with a FSM transition into/out of
      // the current index, and the window component for that transition for
      // the given site in the unit cell.
      std::map<EAIndex, StateComponent> GetWNext(EAIndex i, int n, bool Upper);
      std::map<EAIndex, StateComponent> GetWPrev(EAIndex i, int n, bool Upper);

      // Recurisve methods to get the multiple transtitions for corner elements.
      void GetWNextCorner(std::map<EAIndex, StateComponent>& Result, EAIndex i, int n, bool Upper);
      void GetWPrevCorner(std::map<EAIndex, StateComponent>& Result, EAIndex i, int n, bool Upper);

      // Calculate the unit cell for element (i, j).
      void CalculateE(EAIndex i, EAIndex j);
      void CalculateF(EAIndex i, EAIndex j);

      // Solve the corner element (i, j), using the off-diagonal component CTriK.
      void SolveE(CornerIndex i, CornerIndex j, std::vector<KMatrixPolyType> CTriK);
      void SolveF(CornerIndex i, CornerIndex j, std::vector<KMatrixPolyType> CTriK);

      InfiniteMPO Op;

      // The number of windows in the upper and lower states.
      int NUpper, NLower;
      int Degree;
      double Tol;
      double UnityEpsilon;
      bool NeedFinalMatrix;
      bool EAOptimization;
      bool SubtractEnergy;
      bool PhaseWarnings;
      int Verbose;

      std::map<CornerIndex, LinearWavefunction> PsiUpper, PsiLower;
      std::map<CornerIndex, std::complex<double>> ExpIKUpper, ExpIKLower;
      std::map<std::pair<CornerIndex, CornerIndex>, MatrixOperator> TLeft, TRight;
      std::map<std::pair<CornerIndex, CornerIndex>, bool> TCalculated;
      std::map<EAIndex, std::map<int, StateComponent>> WindowUpper, WindowLower;
      std::map<EAIndex, int> IMax, JMax;
      std::map<std::pair<EAIndex, EAIndex>, std::map<int, std::vector<KMatrixPolyType>>> EMatK, FMatK;
      // If we update the windows, we keep the old E/F matrices as initial guesses for the linear solver.
      std::map<std::pair<EAIndex, EAIndex>, std::vector<KMatrixPolyType>> EMatKOld, FMatKOld;

      int UnitCellSize = 0;
      QuantumNumber QShift;
};

#endif
