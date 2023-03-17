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
#include "mp-algorithms/triangular_mpo_solver_helpers.h"
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

// Shifts the variables in the polynomials in EMatK from (n+Shift) to n.
void
ShiftVariable(std::vector<KMatrixPolyType>& EMatK, int Shift)
{
   // Loop over each element in the vector.
   for (auto& EK : EMatK)
   {
      // Loop over each momentum.
      for (auto& EX : EK)
      {
         // Loop over the terms in the polynomial, starting from the highest degree.
         auto E = EX.second.end();
         while (E != EX.second.begin())
         {
            --E;
            for (int n = 0; n < E->first; ++n)
               EX.second[n] -= double(std::pow(Shift, E->first-n) * Binomial(E->first, n)) * E->second;
         }
      }
   }
}

// Overloaded comparisons for ExtendedInt.
// I haven't defined comparisons with int, since they can convert implcitly to ExtendedInt.
bool
operator==(ExtendedInt I1, ExtendedInt I2)
{
   if (I1.is_inf() || I2.is_inf())
      return I1.is_inf() && I2.is_inf();
   else
      return I1.value() == I2.value();
}

bool
operator<(ExtendedInt I1, ExtendedInt I2)
{
   if (I1.is_inf() || I2.is_inf())
      return I2.is_inf() && !I1.is_inf();
   else
      return I1.value() < I2.value();
}

bool
operator>(ExtendedInt I1, ExtendedInt I2)
{
   if (I1.is_inf() || I2.is_inf())
      return I1.is_inf() && !I2.is_inf();
   else
      return I1.value() > I2.value();
}

bool
operator<=(ExtendedInt I1, ExtendedInt I2)
{
   if (I1.is_inf() || I2.is_inf())
      return I2.is_inf();
   else
      return I1.value() <= I2.value();
}

bool
operator>=(ExtendedInt I1, ExtendedInt I2)
{
   if (I1.is_inf() || I2.is_inf())
      return I1.is_inf();
   else
      return I1.value() >= I2.value();
}

EFMatrix::EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings)
   : Op(Op_), Degree(Settings.Degree), Tol(Settings.Tol),
     UnityEpsilon(Settings.UnityEpsilon), NeedFinalMatrix(Settings.NeedFinalMatrix),
     EAOptimization(Settings.EAOptimization), SubtractEnergy(Settings.SubtractEnergy),
     Verbose(Settings.Verbose)
{
   IMax = 1;
   JMax = 1;

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
EFMatrix::SetPsi(bool i, InfiniteWavefunctionLeft const& Psi, std::complex<double> ExpIK)
{
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

   // Set the momentum factors.
   ExpIKUpper[i] = ExpIK;
   ExpIKLower[i] = ExpIK;

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
EFMatrix::SetDiagTEVsLC(bool i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TRight[std::make_pair(i, i)] = delta_shift(Lambda*Lambda, QShift);
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetPsi(bool i, InfiniteWavefunctionRight const& Psi, std::complex<double> ExpIK)
{
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

   // Set the momentum factors.
   ExpIKUpper[i] = ExpIK;
   ExpIKLower[i] = ExpIK;

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
EFMatrix::SetDiagTEVsRC(bool i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = Lambda*Lambda;
   TRight[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetExpIKUpper(std::complex<double> ExpIK)
{
   ExpIKUpper[true] = ExpIK;
}

void
EFMatrix::SetExpIKLower(std::complex<double> ExpIK)
{
   ExpIKLower[true] = ExpIK;
}

void
EFMatrix::SetWindowUpper(std::vector<LinearWavefunction> const& WindowVec)
{
   CHECK_EQUAL(UnitCellSize, WindowVec.size());

   // We assume that all windows have the same size.
   IMax = WindowVec.front().size();

   int n = 0;
   for (auto const& Window : WindowVec)
   {
      int i = 1;
      for (auto const& B : Window)
         WindowUpper[i++][n] = B;
      ++n;
   }

   // Invalidate elements which depend on this window.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      if (I->first.first > 0)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      if (!I->first.first.is_inf())
         I = FMatK.erase(I);
      else
         ++I;
   }
}

void
EFMatrix::SetWindowLower(std::vector<LinearWavefunction> const& WindowVec)
{
   CHECK_EQUAL(UnitCellSize, WindowVec.size());

   // We assume that all windows have the same size.
   JMax = WindowVec.front().size();

   int n = 0;
   for (auto const& Window : WindowVec)
   {
      int j = 1;
      for (auto const& B : Window)
         WindowLower[j++][n] = B;
      ++n;
   }

   // Invalidate elements which depend on this window.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      if (I->first.second > 0)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      if (!I->first.second.is_inf())
         I = FMatK.erase(I);
      else
         ++I;
   }
}

void
EFMatrix::SetWindowUpper(int i, std::deque<StateComponent> const& BDeque)
{
   PRECONDITION(i > 0)(i <= IMax);
   CHECK_EQUAL(UnitCellSize, BDeque.size());

   int n = 0;
   for (auto const& B : BDeque)
      WindowUpper[i][n++] = B;

   // Invalidate elements which depend on this window.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      if (I->first.first >= i)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      if (I->first.first < i)
         I = FMatK.erase(I);
      else
         ++I;
   }
}

void
EFMatrix::SetWindowLower(int j, std::deque<StateComponent> const& BDeque)
{
   PRECONDITION(j > 0)(j <= JMax);
   CHECK_EQUAL(UnitCellSize, BDeque.size());

   int n = 0;
   for (auto const& B : BDeque)
      WindowLower[j][n++] = B;

   // Invalidate elements which depend on this window element.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      if (I->first.second >= j)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      if (I->first.second < j)
         I = FMatK.erase(I);
      else
         ++I;
   }
}

StateComponent
EFMatrix::GetWUpper(int i, int n)
{
   // If not defined, return null.
   if (WindowUpper.count(i) == 0)
      return StateComponent();

   StateComponent WUpper = WindowUpper[i][numerics::divp(n-i+1,UnitCellSize).rem];
   for (int I = numerics::divp(n-i+1,UnitCellSize).quot; I < 0; ++I)
      WUpper.delta_shift(QShift);
   return WUpper;
}

StateComponent
EFMatrix::GetWLower(int j, int n)
{
   // If not defined, return null;
   if (WindowLower.count(j) == 0)
      return StateComponent();

   StateComponent WLower = WindowLower[j][numerics::divp(n-j+1,UnitCellSize).rem];
   for (int I = numerics::divp(n-j+1,UnitCellSize).quot; I < 0; ++I)
      WLower.delta_shift(QShift);
   return WLower;
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
EFMatrix::GetTLeft(bool i, bool j, bool F)
{
   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   if (F)
      return delta_shift(TLeft[std::make_pair(i, j)], adjoint(QShift));
   else
      return TLeft[std::make_pair(i, j)];
}

MatrixOperator
EFMatrix::GetTRight(bool i, bool j, bool F)
{
   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   if (F)
      return delta_shift(TRight[std::make_pair(i, j)], adjoint(QShift));
   else
      return TRight[std::make_pair(i, j)];
}

void
EFMatrix::CalculateTEVs(bool i, bool j)
{
   ProductMPO StringOp = ProductMPO::make_identity(ExtractLocalBasis(PsiUpper[i]));
   if (Op.is_product())
      StringOp = Op.as_product_mpo();

   // Find the eigenpair if they have an eigenvalue of magnitude 1.
   std::tie(std::ignore, TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)])
      = get_transfer_unit_eigenpair(PsiUpper[i], PsiLower[j], QShift, StringOp, Tol, UnityEpsilon, Verbose-1);

   TRight[std::make_pair(i, j)].delta_shift(QShift);

   // Special normalization for the first transfer matrix.
   if (i == false && j == false)
      Normalize(TRight[std::make_pair(i, j)], TLeft[std::make_pair(i, j)]);

   // Special normalization for the final transfer matrix.
   if (i == true && j == true)
      Normalize(TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)]);

   TCalculated[std::make_pair(i, j)] = true;
}

std::vector<KMatrixPolyType>
EFMatrix::GetE(ExtendedInt i, ExtendedInt j, int n)
{
   // Calculate the element if we haven't already.
   if (EMatK[std::make_pair(i, j)].empty())
      this->CalculateE(i, j);

   if (n == -1)
   {
      // Delta shift and multiply by momentum factor.
      std::vector<KMatrixPolyType> Result = ScalarMultiply(this->MomentumFactor(i, j), delta_shift(EMatK[std::make_pair(i, j)][UnitCellSize-1], QShift));
      // Shift polynomial variable from n+1 to n.
      ShiftVariable(Result, 1);
      return Result;
   }
   else if (n >= 0 && n < UnitCellSize)
      return EMatK[std::make_pair(i, j)][n];
   else
      throw std::runtime_error("EFMatrix::GetE: fatal: n is out of bounds.");
}

std::map<ExtendedInt, StateComponent>
EFMatrix::GetWUpperPrev(ExtendedInt i, int n)
{
   std::map<ExtendedInt, StateComponent> Result;

   if (i > 0)
   {
      if (!i.is_inf())
         Result[i.value()-1] = this->GetWUpper(i.value(), n);
      else
         Result[IMax-1] = this->GetWUpper(IMax, n);
   }

   return Result;
}

std::map<ExtendedInt, StateComponent>
EFMatrix::GetWLowerPrev(ExtendedInt j, int n)
{
   std::map<ExtendedInt, StateComponent> Result;

   if (j > 0)
   {
      if (!j.is_inf())
         Result[j.value()-1] = this->GetWLower(j.value(), n);
      else
         Result[JMax-1] = this->GetWLower(JMax, n);
   }

   return Result;
}

void
EFMatrix::CalculateE(ExtendedInt i, ExtendedInt j)
{
   auto O = Op.as_generic_mpo().begin();
   LinearWavefunction::const_iterator CUpper, CLower;

   bool CornerUpper = i == 0 || i.is_inf();
   if (CornerUpper)
      CUpper = PsiUpper[i.is_inf()].begin();

   bool CornerLower = j == 0 || j.is_inf();
   if (CornerLower)
      CLower = PsiLower[j.is_inf()].begin();

   // Cumulative sum for corner elements.
   std::vector<KMatrixPolyType> CTriK = std::vector<KMatrixPolyType>(O->Basis2().size());

   // Loop over each position in the unit cell.
   for (int n = 0; n < UnitCellSize; ++n)
   {
      EMatK[std::make_pair(i, j)][n] = std::vector<KMatrixPolyType>(O->Basis2().size());

      // Handle contributions with a window on the top and the bottom.
      for (auto const& I : this->GetWUpperPrev(i, n))
         for (auto const& J : this->GetWLowerPrev(j, n))
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(I.second), this->GetE(I.first, J.first, n-1), J.second);

      // Handle contributions with a window on the bottom only.
      if (CornerUpper)
         for (auto const& J : this->GetWLowerPrev(j, n))
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(*CUpper), this->GetE(i, J.first, n-1), J.second);

      // Handle contributions with a window on the top only.
      if (CornerLower)
         for (auto const& I : this->GetWUpperPrev(i, n))
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(I.second), this->GetE(I.first, j, n-1), *CLower);

      // Handle contributions without any windows (corner elements only).
      if (CornerUpper && CornerLower)
      {
         if (n > 0)
            CTriK = contract_from_left(*O, herm(*CUpper), CTriK, *CLower);
         CTriK += EMatK[std::make_pair(i, j)][n];
      }

      ++O;
      if (CornerUpper)
         ++CUpper;
      if (CornerLower)
         ++CLower;
   }

   // Run the linear solver for corner elements.
   if (CornerUpper && CornerLower)
      this->SolveE(i.is_inf(), j.is_inf(), delta_shift(CTriK, QShift));
}

void
EFMatrix::SolveE(bool i, bool j, std::vector<KMatrixPolyType> CTriK)
{
   EMatK[std::make_pair(i, j)][UnitCellSize-1].clear();

   bool FirstElement = i == false && j == false;
   bool FinalElement = i == true && j == true;

   // Initialize the first element of the first E matrix.
   if (FirstElement)
   {
      EMatK[std::make_pair(i, j)][UnitCellSize-1].push_back(KMatrixPolyType());
      EMatK[std::make_pair(i, j)][UnitCellSize-1][0][1.0] = MatrixPolyType(this->GetTLeft(i, j));
   }

   // Run the linear solver.
   SolveMPO_EA_Left(EMatK[std::make_pair(i, j)][UnitCellSize-1], CTriK, PsiUpper[i], PsiLower[j], QShift,
                    Op, this->GetTLeft(i, j), this->GetTRight(i, j), this->MomentumFactor(i, j),
                    Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                    !FirstElement && EAOptimization, Verbose-1);

   EMatK[std::make_pair(i, j)][UnitCellSize-1] = delta_shift(EMatK[std::make_pair(i, j)][UnitCellSize-1], adjoint(QShift));

   // Shift polynomial variable from n-1 to n.
   // TODO: Can this be factored out?
   ShiftVariable(EMatK[std::make_pair(i, j)][UnitCellSize-1], -1);

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
EFMatrix::GetF(ExtendedInt i, ExtendedInt j, int n)
{
   // Calculate the element if we haven't already.
   if (FMatK[std::make_pair(i, j)].empty())
      this->CalculateF(i, j);

   if (n == UnitCellSize)
   {
      // Delta shift and multiply by momentum factor.
      std::vector<KMatrixPolyType> Result = ScalarMultiply(std::conj(this->MomentumFactor(i, j)), delta_shift(FMatK[std::make_pair(i, j)][0], adjoint(QShift)));
      // Shift polynomial variable from n+1 to n.
      ShiftVariable(Result, 1);
      return Result;
   }
   else if (n >= 0 && n < UnitCellSize)
      return FMatK[std::make_pair(i, j)][n];
   else
      throw std::runtime_error("EFMatrix::GetF: fatal: n is out of bounds.");
}

std::map<ExtendedInt, StateComponent>
EFMatrix::GetWUpperNext(ExtendedInt i, int n)
{
   std::map<ExtendedInt, StateComponent> Result;

   if (!i.is_inf())
   {
      if (i.value() < IMax-1)
         Result[i.value()+1] = this->GetWUpper(i.value()+1, n);
      else
         Result[true] = this->GetWUpper(IMax, n);
   }

   return Result;
}

std::map<ExtendedInt, StateComponent>
EFMatrix::GetWLowerNext(ExtendedInt j, int n)
{
   std::map<ExtendedInt, StateComponent> Result;

   if (!j.is_inf())
   {
      if (j.value() < JMax-1)
         Result[j.value()+1] = this->GetWLower(j.value()+1, n);
      else
         Result[true] = this->GetWLower(JMax, n);
   }

   return Result;
}

void
EFMatrix::CalculateF(ExtendedInt i, ExtendedInt j)
{
   auto O = Op.as_generic_mpo().end();
   LinearWavefunction::const_iterator CUpper, CLower;

   bool CornerUpper = i == 0 || i.is_inf();
   if (CornerUpper)
      CUpper = PsiUpper[i.is_inf()].end();

   bool CornerLower = j == 0 || j.is_inf();
   if (CornerLower)
      CLower = PsiLower[j.is_inf()].end();

   // Cumulative sum for corner elements.
   std::vector<KMatrixPolyType> CTriK = std::vector<KMatrixPolyType>((O-1)->Basis1().size());

   // Loop over each position in the unit cell.
   for (int n = UnitCellSize-1; n >= 0; --n)
   {
      --O;
      if (CornerUpper)
         --CUpper;
      if (CornerLower)
         --CLower;

      FMatK[std::make_pair(i, j)][n] = std::vector<KMatrixPolyType>(O->Basis1().size());

      // Handle contributions with a window on the top and the bottom.
      for (auto const& I : this->GetWUpperNext(i, n))
         for (auto const& J : this->GetWLowerNext(j, n))
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), I.second, this->GetF(I.first, J.first, n+1), herm(J.second));

      // Handle contributions with a window on the bottom only.
      if (CornerUpper)
         for (auto const& J : this->GetWLowerNext(j, n))
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), *CUpper, this->GetF(i, J.first, n+1), herm(J.second));

      // Handle contributions with a window on the top only.
      if (CornerLower)
         for (auto const& I : this->GetWUpperNext(i, n))
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), I.second, this->GetF(I.first, j, n+1), herm(*CLower));

      // Handle contributions without any windows (corner elements only).
      if (CornerUpper && CornerLower)
      {
         if (n < UnitCellSize-1)
            CTriK = contract_from_right(herm(*O), *CUpper, CTriK, herm(*CLower));
         CTriK += FMatK[std::make_pair(i, j)][n];
      }
   }

   // Run the linear solver for corner elements.
   if (CornerUpper && CornerLower)
      this->SolveF(i.is_inf(), j.is_inf(), delta_shift(CTriK, adjoint(QShift)));
}

void
EFMatrix::SolveF(bool i, bool j, std::vector<KMatrixPolyType> CTriK)
{
   FMatK[std::make_pair(i, j)][0].clear();

   bool FirstElement = i == true && j == true;
   bool FinalElement = i == false && j == false;

   // Initialize the first element of the first F matrix.
   if (FirstElement)
   {
      FMatK[std::make_pair(i, j)][0].push_back(KMatrixPolyType());
      FMatK[std::make_pair(i, j)][0][0][1.0] = MatrixPolyType(this->GetTRight(i, j, true));
   }

   // Run the linear solver.
   SolveMPO_EA_Right(FMatK[std::make_pair(i, j)][0], CTriK, PsiUpper[i], PsiLower[j], QShift,
                     Op, this->GetTLeft(i, j, true), this->GetTRight(i, j, true), std::conj(this->MomentumFactor(i, j)),
                     Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                     !FirstElement && EAOptimization, Verbose-1);

   FMatK[std::make_pair(i, j)][0] = delta_shift(FMatK[std::make_pair(i, j)][0], QShift);

   // Subtract the contribution due to the energy density: this is to ensure
   // that the expectation value obtained by taking product of E and F matrices
   // corresponds to the value we would obtain by solving the full E or F
   // matrix.
   if (FirstElement && SubtractEnergy)
   {
      // The first contribution can be extracted from the expectation value of the F matrix.
      std::complex<double> RightEnergy = inner_prod(this->GetTLeft(i, j), FMatK[std::make_pair(i, j)][0].front()[1.0].coefficient(1));

      if (Verbose > 0)
         std::cerr << "EFMatrix: subtracting contribution from right energy = " << RightEnergy << std::endl;

      // The second contribution is the "bond energy", which is the
      // contribution to the energy from terms which cross the unit cell
      // boundary: the most reliable way to calculate this is by calculating
      // the E matrix for the same wavefunctions and taking the product of
      // the E and F matrices.
      std::vector<KMatrixPolyType> EMatKRight(1, KMatrixPolyType());
      EMatKRight[0][1.0] = MatrixPolyType(this->GetTLeft(i, j));

      SolveMPO_EA_Left(EMatKRight, std::vector<KMatrixPolyType>(), PsiUpper[i], PsiLower[j], QShift,
                       Op, this->GetTLeft(i, j), this->GetTRight(i, j), std::conj(this->MomentumFactor(i, j)),
                       Degree, Tol, UnityEpsilon, true, false, Verbose-1);

      std::complex<double> BondEnergy = 0.0;
      for (int I = 0; I < EMatKRight.size(); ++I)
         BondEnergy += inner_prod(EMatKRight[I][1.0].coefficient(0), FMatK[std::make_pair(i, j)][0][I][1.0].coefficient(0));

      if (Verbose > 0)
         std::cerr << "EFMatrix: subtracting contribution from bond energy = " << BondEnergy << std::endl;

      // Subtract the right and bond energies from the final element from the first F matrix.
      FMatK[std::make_pair(i, j)][0].front() -= (RightEnergy + BondEnergy) * FMatK[std::make_pair(i, j)][0].back();
   }

   // Shift polynomial variable from n-1 to n.
   // TODO: Can this be factored out?
   ShiftVariable(FMatK[std::make_pair(i, j)][0], -1);

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

// TODO: Could these be implemented in a better way?
StateComponent
EFMatrix::GetESC(ExtendedInt i, ExtendedInt j, int n)
{
   std::vector<KMatrixPolyType> E = this->GetE(i, j, n);
   MatrixOperator Tmp = E[0][1.0].coefficient(0);
   OperatorComponent O = Op.as_generic_mpo()[numerics::divp(n, UnitCellSize).rem];
   StateComponent Result(O.Basis2(), Tmp.Basis1(), Tmp.Basis2());

   for (int I = 0; I < Result.size(); ++I)
      Result[I] = E[I][1.0].coefficient(0);

   return Result;
}

StateComponent
EFMatrix::GetFSC(ExtendedInt i, ExtendedInt j, int n)
{
   std::vector<KMatrixPolyType> F = this->GetF(i, j, n);
   MatrixOperator Tmp = F[0][1.0].coefficient(0);
   OperatorComponent O = Op.as_generic_mpo()[numerics::divp(n, UnitCellSize).rem];
   StateComponent Result(O.Basis1(), Tmp.Basis1(), Tmp.Basis2());

   for (int I = 0; I < Result.size(); ++I)
      Result[I] = F[I][1.0].coefficient(0);

   return Result;
}

std::deque<StateComponent>
EFMatrix::GetHEff(int i)
{
   std::deque<StateComponent> Result;

   auto CLeft = PsiLower[false].begin();
   auto CRight = PsiLower[true].begin();
   auto O = Op.as_generic_mpo().begin();

   ExtendedInt iNext;
   if (i == IMax-1)
      iNext = true;
   else
      iNext = i+1;

   // Loop over each position in the unit cell.
   for (int n = 0; n < UnitCellSize; ++n)
   {
      // Add each contribution due to each component on the bottom.
      Result.push_back(operator_prod_inner(*O, this->GetESC(i, 0, n-1), *CLeft, herm(this->GetFSC(iNext, 0, n+1))));

      for (int j = 0; j < JMax; ++j)
      {
         ExtendedInt jNext;
         if (j == JMax-1)
            jNext = true;
         else
            jNext = j+1;

         Result.back() += operator_prod_inner(*O, this->GetESC(i, j, n-1), this->GetWLower(j+1, n), herm(this->GetFSC(iNext, jNext, n+1)));
      }

      Result.back() += operator_prod_inner(*O, this->GetESC(i, true, n-1), *CRight, herm(this->GetFSC(iNext, true, n+1)));

      ++CLeft, ++CRight, ++O;
   }

   // Rotate the result deque so the components are in window order, not unit cell order.
   for (int I = 0; I < i; ++I)
   {
      Result.push_back(delta_shift(Result.front(), adjoint(QShift)));
      Result.pop_front();
   }

   return Result;
}
