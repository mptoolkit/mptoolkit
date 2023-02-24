// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ef-matrix.cpp
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

#include "ef-matrix.h"
#include "tensor/tensor_eigen.h"

// The tolerance of the trace of the left/right generalised transfer
// eigenvectors for fixing their relative phases.
double const TraceTol = 1e-8;

// Normalize two TEVs such that the sum of singular values of Rho is one and
// the trace of Ident has phase zero: this important if we cannot get Rho and
// Ident from the lambda matrix (e.g. if we are using a string operator or if
// PsiUpper and PsiLower are different). This is to ensure the final
// expectation value has the correct magnitude and to remove any spurious phase
// due to the eigensolver.
void
Normalize(MatrixOperator& Rho, MatrixOperator& Ident)
{
   if (Rho.is_null())
      throw std::runtime_error("EFMatrix::CalculateTEVs: fatal: leading transfer matrix eigenvalue is below tolerance for diagonal element");

   // Normalize by setting the sum of singular values of Rho to be 1.
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SingularValueDecomposition(Rho, U, D, Vh);
   Rho *= 1.0 / trace(D);
   Ident *= 1.0 / inner_prod(Rho, Ident);

   // Fix the phases of Rho and Ident by setting the phase of the trace of Ident to be zero.
   std::complex<double> Trace = trace(Ident);
   if (std::abs(Trace) > TraceTol)
   {
      std::complex<double> ConjPhase = std::conj(Trace) / std::abs(Trace);
      Ident *= ConjPhase;
      Rho *= ConjPhase;
   }
   else
      std::cerr << "EFMatrix: warning: the trace of Ident is below threshold,"
                   " so the results will have a spurious phase contribution." << std::endl;
}

EFMatrix::EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings)
   : Op(Op_), Degree(Settings.Degree), Tol(Settings.Tol),
     UnityEpsilon(Settings.UnityEpsilon), NeedFinalMatrix(Settings.NeedFinalMatrix),
     EAOptimization(Settings.EAOptimization), Verbose(Settings.Verbose)
{
   ExpIKUpper[0] = 1.0;
   ExpIKLower[0] = 1.0;
   IMax = -1;
   JMax = -1;

   this->CheckOperator();
}

void
EFMatrix::CheckOperator()
{
   if (Op.is_product())
   {
      // We only handle string operators at the moment.
      if (!Op.as_product_mpo().is_string())
         throw std::runtime_error("EFMatrix: fatal: Op cannot be a non-string ProductMPO.");
      // We only handle scalar operators at the moment.
      if (!Op.as_product_mpo().is_scalar())
         throw std::runtime_error("EFMatrix: fatal: Op cannot be a non-scalar ProductMPO.");
   }
   else if (Op.is_triangular())
   {
      // TODO: We could possibly relax this.
      if (!classify(Op.as_basic_triangular_mpo()(0,0), UnityEpsilon).is_identity())
         throw std::runtime_error("EFMatrix: fatal: first component of Op must be the identity operator.");
      int Dim = Op.as_basic_triangular_mpo().Basis1().size();
      if (!classify(Op.as_basic_triangular_mpo()(Dim-1,Dim-1), UnityEpsilon).is_identity())
         throw std::runtime_error("EFMatrix: fatal: final component of Op must be the identity operator.");
   }
   else
      throw std::runtime_error("EFMatrix: fatal: Op must be a BasicTriangularMPO or a ProductMPO.");
}

void
EFMatrix::SetPsi(int i, InfiniteWavefunctionLeft const& Psi)
{
   PRECONDITION(i >= 0);

   // If we haven't set QShift, set it now, otherwise, check that it is the same.
   if (QShift.is_null())
      QShift = Psi.qshift();
   else
      CHECK_EQUAL(QShift, Psi.qshift());

   // If we haven't set the unit cell size, set it now, otherwise, check that it is the same.
   if (UnitCellSize == 0)
      UnitCellSize = Psi.size();
   else
      CHECK_EQUAL(UnitCellSize, Psi.size());

   // Update the maximum indices.
   if (i > IMax)
      IMax = i;
   if (i > JMax)
      JMax = i;

   // Extract the wavefunction.
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(PsiLinear, Lambda) = get_left_canonical(Psi);

   PsiUpper[i] = PsiLinear;
   PsiLower[i] = PsiLinear;

   // Set the transfer matrix eigenvectors if we have a triangular operator or the identity.
   if (!Op.is_product() || Op.as_product_mpo().is_identity())
      this->SetDiagTEVsLC(i, Lambda);
}

void
EFMatrix::SetDiagTEVsLC(int i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TRight[std::make_pair(i, i)] = delta_shift(Lambda*Lambda, QShift);
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetPsi(int i, InfiniteWavefunctionRight const& Psi)
{
   PRECONDITION(i >= 0);

   // If we haven't set QShift, set it now, otherwise, check that it is the same.
   if (QShift.is_null())
      QShift = Psi.qshift();
   else
      CHECK_EQUAL(QShift, Psi.qshift());

   // If we haven't set the unit cell size, set it now, otherwise, check that it is the same.
   if (UnitCellSize == 0)
      UnitCellSize = Psi.size();
   else
      CHECK_EQUAL(UnitCellSize, Psi.size());

   // Update the maximum indices.
   if (i > IMax)
      IMax = i;
   if (i > JMax)
      JMax = i;

   // Extract the wavefunction.
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(Lambda, PsiLinear) = get_right_canonical(Psi);

   PsiUpper[i] = PsiLinear;
   PsiLower[i] = PsiLinear;

   // Set the transfer matrix eigenvectors if we have a triangular operator or the identity.
   if (!Op.is_product() || Op.as_product_mpo().is_identity())
      this->SetDiagTEVsRC(i, Lambda);
}

void
EFMatrix::SetDiagTEVsRC(int i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = Lambda*Lambda;
   TRight[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK)
{
   PRECONDITION(i > 0);
   CHECK_EQUAL(UnitCellSize, BDeque.size());

   int n = 0;
   for (auto const& B : BDeque)
   {
      WindowUpper[n][i] = B;
      ++n;
   }

   // We must have already set i-1.
   // FIXME: Is it possible to do this in another way?
   CHECK(ExpIKUpper.count(i-1) == 1);
   ExpIKUpper[i] = ExpIKUpper[i-1] * ExpIK;

   // TODO: Invalidation.
}

void
EFMatrix::SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK)
{
   PRECONDITION(j > 0);
   CHECK_EQUAL(UnitCellSize, BDeque.size());

   int n = 0;
   for (auto const& B : BDeque)
   {
      WindowLower[n][j] = B;
      ++n;
   }

   // We must have already set j-1.
   // FIXME: Is it possible to do this in another way?
   CHECK(ExpIKLower.count(j-1) == 1);
   ExpIKLower[j] = ExpIKLower[j-1] * ExpIK;
}

void
EFMatrix::SetOp(InfiniteMPO Op_, int Degree_)
{
   // If either the new or old Op is a product operator except for the
   // identity, we need to recalculate the transfer matrix eigenvalues.
   if ((Op.is_product() && !Op.as_product_mpo().is_identity())
      || (Op_.is_product() && !Op_.as_product_mpo().is_identity()))
   {
      TCalculated.clear();
      TLeft.clear();
      TRight.clear();
   }
   // TODO: If the new operator is triangular or the identity, and we have the
   // Lambdas for PsiUpper/Lower, we should be able to get the EVs without
   // using the solver.

   // Clear the matrix elements.
   EMatK.clear();
   FMatK.clear();

   Op = Op_;
   Degree = Degree_;

   // TODO: Check Op size.

   this->CheckOperator();
}

MatrixOperator
EFMatrix::GetTLeft(int i, int j, bool F)
{
   // Return null if out of bounds.
   if (i < 0 || j < 0)
      return MatrixOperator();

   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   if (F)
      return delta_shift(TLeft[std::make_pair(i, j)], adjoint(QShift));
   else
      return TLeft[std::make_pair(i, j)];
}

MatrixOperator
EFMatrix::GetTRight(int i, int j, bool F)
{
   // Return null if out of bounds.
   if (i < 0 || j < 0)
      return MatrixOperator();

   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   if (F)
      return delta_shift(TRight[std::make_pair(i, j)], adjoint(QShift));
   else
      return TRight[std::make_pair(i, j)];
}

void
EFMatrix::CalculateTEVs(int i, int j)
{
   ProductMPO StringOp = ProductMPO::make_identity(ExtractLocalBasis(PsiUpper[i]));
   if (Op.is_product())
      StringOp = Op.as_product_mpo();

   // Find the eigenpair if they have an eigenvalue of magnitude 1.
   std::tie(std::ignore, TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)])
      = get_transfer_unit_eigenpair(PsiUpper[i], PsiLower[j], QShift, StringOp, Tol, UnityEpsilon, Verbose);

   TRight[std::make_pair(i, j)].delta_shift(QShift);

   // Special normalization for the first transfer matrix.
   if (i == 0 && j == 0)
      Normalize(TRight[std::make_pair(i, j)], TLeft[std::make_pair(i, j)]);

   // Special normalization for the final transfer matrix.
   if (i == IMax && j == JMax)
      Normalize(TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)]);

   TCalculated[std::make_pair(i, j)] = true;
}

std::vector<KMatrixPolyType>
EFMatrix::GetE(int i, int j, int n)
{
   // Return null if out of bounds.
   if (i < 0 || j < 0)
      return std::vector<KMatrixPolyType>();

   // Calculate the element if we haven't already.
   if (EMatK[std::make_pair(i, j)].empty())
      this->CalculateE(i, j);

   if (n == -1)
      return ScalarMultiply(this->MomentumFactor(i, j), delta_shift(EMatK[std::make_pair(i, j)][UnitCellSize-1], QShift));
   else if (n == UnitCellSize)
      return ScalarMultiply(std::conj(this->MomentumFactor(i, j)), delta_shift(EMatK[std::make_pair(i, j)][0], adjoint(QShift)));
   else if (n >= 0 && n < UnitCellSize)
      return EMatK[std::make_pair(i, j)][n];
   else
      throw std::runtime_error("EFMatrix::GetE: fatal: n is out of bounds.");
}

void
EFMatrix::CalculateE(int i, int j)
{
   auto CUpper = PsiUpper[i].begin();
   auto CLower = PsiLower[j].begin();
   auto O = Op.as_generic_mpo().begin();

   // Cumulative sum for corner elements.
   std::vector<KMatrixPolyType> CTriK = std::vector<KMatrixPolyType>(O->Basis2().size());

   // Loop over each position in the unit cell.
   for (int n = 0; n < UnitCellSize; ++n)
   {
      EMatK[std::make_pair(i, j)][n] = std::vector<KMatrixPolyType>(O->Basis2().size());

      // Handle contributions with a window on the top and the bottom.
      if (WindowUpper[numerics::divp(n-i+1,UnitCellSize).rem].count(i) == 1)
         if (WindowLower[numerics::divp(n-j+1,UnitCellSize).rem].count(j) == 1)
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(WindowUpper[numerics::divp(n-i+1,UnitCellSize).rem][i]), this->GetE(i-1, j-1, n-1), WindowLower[numerics::divp(n-j+1,UnitCellSize).rem][j]);

      // Handle contributions with a window on the bottom only.
      if (i == 0 || i == IMax)
         if (WindowLower[numerics::divp(n-j+1,UnitCellSize).rem].count(j) == 1)
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(*CUpper), this->GetE(i, j-1, n-1), WindowLower[numerics::divp(n-j+1,UnitCellSize).rem][j]);

      // Handle contributions with a window on the top only.
      if (j == 0 || j == JMax)
         if (WindowUpper[numerics::divp(n-i+1,UnitCellSize).rem].count(i) == 1)
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(WindowUpper[numerics::divp(n-i+1,UnitCellSize).rem][i]), this->GetE(i-1, j, n-1), *CLower);

      // Handle contributions without any windows (corner elements only).
      if ((i == 0 || i == IMax) && (j == 0 || j == JMax))
      {
         if (n > 0)
            CTriK = contract_from_left(*O, herm(*CUpper), CTriK, *CLower);
         CTriK += EMatK[std::make_pair(i, j)][n];
      }

      ++CUpper, ++CLower, ++O;
   }

   // Run the linear solver for corner elements.
   if ((i == 0 || i == IMax) && (j == 0 || j == JMax))
      this->SolveE(i, j, delta_shift(CTriK, QShift));
}

void
EFMatrix::SolveE(int i, int j, std::vector<KMatrixPolyType> CTriK)
{
   EMatK[std::make_pair(i, j)][UnitCellSize-1].clear();

   // Initialize the first element of the first E matrix.
   if (i == 0 && j == 0)
   {
      EMatK[std::make_pair(i, j)][UnitCellSize-1].push_back(KMatrixPolyType());
      EMatK[std::make_pair(i, j)][UnitCellSize-1][0][1.0] = MatrixPolyType(this->GetTLeft(i, j));
   }

   bool FinalElement = i == IMax && j == JMax;

   // Run the linear solver.
   SolveMPO_EA_Left(EMatK[std::make_pair(i, j)][UnitCellSize-1], CTriK, PsiUpper[i], PsiLower[j], QShift,
                    Op, this->GetTLeft(i, j), this->GetTRight(i, j), this->MomentumFactor(i, j),
                    Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                    FinalElement && EAOptimization, Verbose);

   EMatK[std::make_pair(i, j)][UnitCellSize-1] = delta_shift(EMatK[std::make_pair(i, j)][UnitCellSize-1], adjoint(QShift));

   // Calculate the elements in the rest of the unit cell.
   if (!FinalElement || NeedFinalMatrix)
   {
      auto CUpper = PsiUpper[i].begin();
      auto CLower = PsiLower[j].begin();
      auto O = Op.as_generic_mpo().begin();

      for (int n = 0; n < UnitCellSize-1; ++n)
      {
         EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(*CUpper), this->GetE(i, j, n-1), *CLower);
         ++CUpper, ++CLower, ++O;
      }
   }
}

std::vector<KMatrixPolyType>
EFMatrix::GetF(int i, int j, int n)
{
   // Return null if out of bounds.
   if (i > IMax || j > JMax)
      return std::vector<KMatrixPolyType>();

   // Calculate the element if we haven't already.
   if (FMatK[std::make_pair(i, j)].empty())
      this->CalculateF(i, j);

   if (n == -1)
      return ScalarMultiply(this->MomentumFactor(i, j), delta_shift(FMatK[std::make_pair(i, j)][UnitCellSize-1], QShift));
   else if (n == UnitCellSize)
      return ScalarMultiply(std::conj(this->MomentumFactor(i, j)), delta_shift(FMatK[std::make_pair(i, j)][0], adjoint(QShift)));
   else if (n >= 0 && n < UnitCellSize)
      return FMatK[std::make_pair(i, j)][n];
   else
      throw std::runtime_error("EFMatrix::GetF: fatal: n is out of bounds.");
}

void
EFMatrix::CalculateF(int i, int j)
{
   auto CUpper = PsiUpper[i].end();
   auto CLower = PsiLower[j].end();
   auto O = Op.as_generic_mpo().end();

   // Cumulative sum for corner elements.
   std::vector<KMatrixPolyType> CTriK = std::vector<KMatrixPolyType>((O-1)->Basis1().size());

   // Loop over each position in the unit cell.
   for (int n = UnitCellSize-1; n >= 0; --n)
   {
      --CUpper, --CLower, --O;

      FMatK[std::make_pair(i, j)][n] = std::vector<KMatrixPolyType>(O->Basis1().size());

      // Handle contributions with a window on the top and the bottom.
      if (WindowUpper[numerics::divp(n-i,UnitCellSize).rem].count(i+1) == 1)
         if (WindowLower[numerics::divp(n-j,UnitCellSize).rem].count(j+1) == 1)
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), WindowUpper[numerics::divp(n-i,UnitCellSize).rem][i+1], this->GetF(i+1, j+1, n+1), herm(WindowLower[numerics::divp(n-j,UnitCellSize).rem][j+1]));

      // Handle contributions with a window on the bottom only.
      if (i == 0 || i == IMax)
         if (WindowLower[numerics::divp(n-j,UnitCellSize).rem].count(j+1) == 1)
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), *CUpper, this->GetF(i, j+1, n+1), herm(WindowLower[numerics::divp(n-j,UnitCellSize).rem][j+1]));

      // Handle contributions with a window on the top only.
      if (j == 0 || j == JMax)
         if (WindowUpper[numerics::divp(n-i,UnitCellSize).rem].count(i+1) == 1)
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), WindowUpper[numerics::divp(n-i,UnitCellSize).rem][i+1], this->GetF(i+1, j, n+1), herm(*CLower));

      // Handle contributions without any windows (corner elements only).
      if ((i == 0 || i == IMax) && (j == 0 || j == JMax))
      {
         if (n < UnitCellSize-1)
            CTriK = contract_from_right(herm(*O), *CUpper, CTriK, herm(*CLower));
         CTriK += FMatK[std::make_pair(i, j)][n];
      }
   }

   // Run the linear solver for corner elements.
   if ((i == 0 || i == IMax) && (j == 0 || j == JMax))
      this->SolveF(i, j, delta_shift(CTriK, adjoint(QShift)));
}

void
EFMatrix::SolveF(int i, int j, std::vector<KMatrixPolyType> CTriK)
{
   FMatK[std::make_pair(i, j)][0].clear();

   // Initialize the first element of the first F matrix.
   if (i == IMax && j == JMax)
   {
      FMatK[std::make_pair(i, j)][0].push_back(KMatrixPolyType());
      FMatK[std::make_pair(i, j)][0][0][1.0] = MatrixPolyType(this->GetTRight(i, j, true));
   }

   bool FinalElement = i == 0 && j == 0;

   // Run the linear solver.
   SolveMPO_EA_Right(FMatK[std::make_pair(i, j)][0], CTriK, PsiUpper[i], PsiLower[j], QShift,
                     Op, this->GetTLeft(i, j, true), this->GetTRight(i, j, true), std::conj(this->MomentumFactor(i, j)),
                     Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                     FinalElement && EAOptimization, Verbose);

   FMatK[std::make_pair(i, j)][0] = delta_shift(FMatK[std::make_pair(i, j)][0], QShift);

   // Calculate the elements in the rest of the unit cell.
   if (!FinalElement || NeedFinalMatrix)
   {
      auto CUpper = PsiUpper[i].end();
      auto CLower = PsiLower[j].end();
      auto O = Op.as_generic_mpo().end();

      for (int n = UnitCellSize-1; n > 0; --n)
      {
         --CUpper, --CLower, --O;
         FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), *CUpper, this->GetF(i, j, n+1), herm(*CLower));
      }
   }
}
