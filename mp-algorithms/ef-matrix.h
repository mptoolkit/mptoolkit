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
   int Degree = 0;
   double Tol = DefaultTol;
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   bool NeedFinalMatrix = true;
   bool EAOptimization = false;
   bool SubtractEnergy = false;
   int Verbose = 0;
};

// A class to represent an extended integer, which can be a finite integer or infinity.
class ExtendedInt
{
   public:
      ExtendedInt() {}

      ExtendedInt(int Value_)
         : Value(Value_), Inf(false) {}

      ExtendedInt(bool Inf_)
         : Value(0), Inf(Inf_) {}

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
      EFMatrix() {}

      EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings);

      // Set boundary wavefunctions: at the moment, we assume that the boundary
      // wavefunctions for the bra and ket are the same, as this allows us to
      // obtain the TEVs from the lambda matrix.
      void SetPsi(bool i, InfiniteWavefunctionLeft const& Psi, std::complex<double> ExpIK = 1.0);
      void SetPsi(bool i, InfiniteWavefunctionRight const& Psi, std::complex<double> ExpIK = 1.0);

      // Set the momentum factor.
      void SetExpIKUpper(std::complex<double> ExpIK);
      void SetExpIKLower(std::complex<double> ExpIK);

      // Set the windows using a vector of LinearWavefunctions.
      void SetWindowUpper(std::vector<LinearWavefunction> const& WindowVec);
      void SetWindowLower(std::vector<LinearWavefunction> const& WindowVec);

      // Set a window element using a deque of single-site windows for each
      // position in the unit cell (in window order). (0 < i < IMax)
      void SetWindowUpper(int i, std::deque<StateComponent> const& BDeque);
      void SetWindowLower(int j, std::deque<StateComponent> const& BDeque);

      // Update the operator: invalidates the calculated E/F matrix elements and
      // (for string operators) the TEVs.
      void SetOp(InfiniteMPO Op_, int Degree_ = 0);

      // Returns the left and right (unit) TEVs, calculating them if they have
      // not been calculated yet. If the spectral radius is < 1, then
      // TLeft/Right will be null.
      MatrixOperator GetTLeft(bool i, bool j, bool F = false);
      MatrixOperator GetTRight(bool i, bool j, bool F = false);

      // Ident is the TEV in the same direction as the E/F matrix.
      MatrixOperator GetIdent(bool i, bool j, bool F = false)
         { return F ? this->GetTRight(i, j, true) : this->GetTLeft(i, j); }
      // Rho is the TEV in the opposite direction as the E/F matrix
      // (calculating the inner product of an E/F matrix element with Rho gives
      // the expectation value).
      MatrixOperator GetRho(bool i, bool j, bool F = false)
         { return F ? this->GetTLeft(i, j, true) : this->GetTRight(i, j); }

      // Get the E/F matrix where the ket and bra wavefunctions go up to the
      // boundaries i and j respectively.
      std::vector<KMatrixPolyType> GetE(ExtendedInt i, ExtendedInt j, int n);
      std::vector<KMatrixPolyType> GetE(ExtendedInt i, ExtendedInt j) { return this->GetE(i, j, -1); }
      std::vector<KMatrixPolyType> GetF(ExtendedInt i, ExtendedInt j, int n);
      std::vector<KMatrixPolyType> GetF(ExtendedInt i, ExtendedInt j) { return this->GetF(i, j, UnitCellSize); }

      std::vector<KMatrixPolyType> GetElement(ExtendedInt i, ExtendedInt j, bool F = false)
         { return F ? this->GetF(i, j) : this->GetE(i, j); }

      // Get the constant, zero-momentum part of this element as a StateComponent.
      StateComponent GetESC(ExtendedInt i, ExtendedInt j, int n);
      StateComponent GetFSC(ExtendedInt i, ExtendedInt j, int n);

      // Calculate action of the effective Hamiltonian on the window w.r.t.
      // site i for each position in the unit cell (in window order).
      std::deque<StateComponent> GetHEff(int i = 0);

   private:
      void CheckOperator();

      // Set the TEVs when PsiUpper and PsiLower are the same and the lambda
      // matrix for the left/right canonical form is known.
      void SetDiagTEVsLC(bool i, RealDiagonalOperator Lambda);
      void SetDiagTEVsRC(bool i, RealDiagonalOperator Lambda);

      // Solve for the TEVs with an eigenvalue of magnitude 1: if the
      // spectral radius is less than one, set the TEVs to null.
      void CalculateTEVs(bool i, bool j);

      // Return the momentum factor corresponding to the element (i, j).
      std::complex<double> MomentumFactor(ExtendedInt i, ExtendedInt j)
         { return ExpIKUpper[i.is_inf()] * std::conj(ExpIKLower[j.is_inf()]); }

      // Return the StateComponent for the window at position n in the unit
      // cell for the index i/j.
      StateComponent GetWUpper(int i, int n);
      StateComponent GetWLower(int j, int n);

      // Get a map of all possible indices with a FSM transition into/out of
      // the current index, and the window component for that transition for
      // the given site in the unit cell.
      std::map<ExtendedInt, StateComponent> GetWUpperPrev(ExtendedInt i, int n);
      std::map<ExtendedInt, StateComponent> GetWLowerPrev(ExtendedInt j, int n);
      std::map<ExtendedInt, StateComponent> GetWUpperNext(ExtendedInt i, int n);
      std::map<ExtendedInt, StateComponent> GetWLowerNext(ExtendedInt j, int n);

      // Calculate the unit cell for element (i, j).
      void CalculateE(ExtendedInt i, ExtendedInt j);
      void CalculateF(ExtendedInt i, ExtendedInt j);

      // Solve the corner element (i, j), using the off-diagonal component CTriK.
      void SolveE(bool i, bool j, std::vector<KMatrixPolyType> CTriK);
      void SolveF(bool i, bool j, std::vector<KMatrixPolyType> CTriK);

      std::map<bool, LinearWavefunction> PsiUpper, PsiLower;
      std::map<ExtendedInt, std::map<int, StateComponent>> WindowUpper, WindowLower;
      std::map<bool, std::complex<double>> ExpIKUpper, ExpIKLower;
      std::map<std::pair<bool, bool>, MatrixOperator> TLeft, TRight;
      std::map<std::pair<bool, bool>, bool> TCalculated;
      std::map<std::pair<ExtendedInt, ExtendedInt>, std::map<int, std::vector<KMatrixPolyType>>> EMatK, FMatK;

      int IMax, JMax;
      int UnitCellSize = 0;
      QuantumNumber QShift;

      InfiniteMPO Op;

      int Degree;
      double Tol;
      double UnityEpsilon;
      bool NeedFinalMatrix;
      bool EAOptimization;
      bool SubtractEnergy;
      int Verbose;
};

#endif
