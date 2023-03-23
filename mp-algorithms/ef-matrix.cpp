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

// Returns true if the input is zero or infinity.
bool
IsCorner(ExtendedInt Input)
{
   return Input == 0 || Input.is_inf();
}

// Returns true if each element of the input vector is zero or infinity.
bool
IsCorner(std::vector<ExtendedInt> Input)
{
   bool Result = true;
   for (auto const& I : Input)
      Result = Result && IsCorner(I);
   return Result;
}

// Convert a vector of ExtendedInts to a vector of ZeroInf objects for whether
// the value is finite or infinite. (NOTE: Nonzero finite values map to "zero",
// so this function maps to the corner index prior to the current index if it
// is not a corner index.)
std::vector<ZeroInf>
ToCorner(std::vector<ExtendedInt> Input)
{
   std::vector<ZeroInf> Result;
   for (auto const& I : Input)
      Result.push_back(ZeroInf(I.is_inf()));
   return Result;
}

EFMatrix::EFMatrix(InfiniteMPO Op_, EFMatrixSettings Settings)
   : Op(Op_), NUpper(Settings.NUpper), NLower(Settings.NLower),
     Degree(Settings.Degree), Tol(Settings.Tol),
     UnityEpsilon(Settings.UnityEpsilon), NeedFinalMatrix(Settings.NeedFinalMatrix),
     EAOptimization(Settings.EAOptimization), SubtractEnergy(Settings.SubtractEnergy),
     Verbose(Settings.Verbose)
{
   // Set default values of IMax and JMax.
   IMax[EAIndex(NUpper, 1)] = 1;
   JMax[EAIndex(NLower, 1)] = 1;

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
EFMatrix::SetPsi(CornerIndex i, InfiniteWavefunctionLeft const& Psi, std::complex<double> ExpIK)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(i.size(), NLower);

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
   this->SetExpIKUpper(i, ExpIK);
   this->SetExpIKLower(i, ExpIK);

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
EFMatrix::SetDiagTEVsLC(CornerIndex i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TRight[std::make_pair(i, i)] = delta_shift(Lambda*Lambda, QShift);
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetPsi(CornerIndex i, InfiniteWavefunctionRight const& Psi, std::complex<double> ExpIK)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(i.size(), NLower);

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
   this->SetExpIKUpper(i, ExpIK);
   this->SetExpIKLower(i, ExpIK);

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
EFMatrix::SetDiagTEVsRC(CornerIndex i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = Lambda*Lambda;
   TRight[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetExpIKUpper(CornerIndex i, std::complex<double> ExpIK)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   ExpIKUpper[i] = ExpIK;
}

void
EFMatrix::SetExpIKLower(CornerIndex j, std::complex<double> ExpIK)
{
   PRECONDITION_EQUAL(j.size(), NLower);
   ExpIKLower[j] = ExpIK;
}

void
EFMatrix::SetWindowUpper(EAIndex i, std::vector<LinearWavefunction> const& WindowVec)
{
   PRECONDITION_EQUAL(UnitCellSize, WindowVec.size());
   PRECONDITION_EQUAL(i.size(), NUpper);
   for (auto const& I : i)
      PRECONDITION(IsCorner(I) || I == 1);

   // We assume that all windows have the same size.
   IMax[i] = WindowVec.front().size();

   int n = 0;
   for (auto const& Window : WindowVec)
   {
      int Site = 1;
      for (auto const& B : Window)
      {
         // Calculate the index for this element.
         EAIndex Index = i;
         for (auto& I : Index)
            if (I == 1)
               I = Site;

         WindowUpper[Index][n] = B;
         ++Site;
      }
      ++n;
   }

   // Invalidate elements which depend on this window.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      bool Erase = true;
      auto iIter = i.begin();
      auto IndexIter = I->first.first.begin();
      while (iIter != i.end())
      {
         if (*iIter == 1 && *IndexIter == 0)
         {
            Erase = false;
            break;
         }
         ++iIter, ++IndexIter;
      }

      if (Erase)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      bool Erase = true;
      auto iIter = i.begin();
      auto IndexIter = I->first.first.begin();
      while (iIter != i.end())
      {
         if (*iIter == 1 && IndexIter->is_inf())
         {
            Erase = false;
            break;
         }
         ++iIter, ++IndexIter;
      }

      if (Erase)
         I = FMatK.erase(I);
      else
         ++I;
   }
}

void
EFMatrix::SetWindowLower(EAIndex j, std::vector<LinearWavefunction> const& WindowVec)
{
   PRECONDITION_EQUAL(UnitCellSize, WindowVec.size());
   PRECONDITION_EQUAL(j.size(), NLower);
   for (auto const& I : j)
      PRECONDITION(IsCorner(I) || I == 1);

   // We assume that all windows have the same size.
   JMax[j] = WindowVec.front().size();

   int n = 0;
   for (auto const& Window : WindowVec)
   {
      int Site = 1;
      for (auto const& B : Window)
      {
         // Calculate the index for this element.
         EAIndex Index = j;
         for (auto& I : Index)
            if (I == 1)
               I = Site;

         WindowLower[Index][n] = B;
         ++Site;
      }
      ++n;
   }

   // Invalidate elements which depend on this window.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      bool Erase = true;
      auto jIter = j.begin();
      auto IndexIter = I->first.second.begin();
      while (jIter != j.end())
      {
         if (*jIter == 1 && *IndexIter == 0)
         {
            Erase = false;
            break;
         }
         ++jIter, ++IndexIter;
      }

      if (Erase)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      bool Erase = true;
      auto jIter = j.begin();
      auto IndexIter = I->first.second.begin();
      while (jIter != j.end())
      {
         if (*jIter == 1 && IndexIter->is_inf())
         {
            Erase = false;
            break;
         }
         ++jIter, ++IndexIter;
      }

      if (Erase)
         I = FMatK.erase(I);
      else
         ++I;
   }
}

void
EFMatrix::SetWUpper(int i, std::deque<StateComponent> const& BDeque)
{
   PRECONDITION(i > 0)(i <= IMax[EAIndex(NUpper, 1)]);
   PRECONDITION_EQUAL(UnitCellSize, BDeque.size());

   EAIndex Index(NUpper, i);

   int n = 0;
   for (auto const& B : BDeque)
      WindowUpper[Index][n++] = B;

   // Invalidate elements which depend on this window element.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      bool Erase = true;
      for (auto const& J : I->first.first)
      {
         if (J < i)
         {
            Erase = false;
            break;
         }
      }

      if (Erase)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      bool Erase = true;
      for (auto const& J : I->first.first)
      {
         if (J >= i)
         {
            Erase = false;
            break;
         }
      }

      if (Erase)
         I = FMatK.erase(I);
      else
         ++I;
   }
}

void
EFMatrix::SetWLower(int j, std::deque<StateComponent> const& BDeque)
{
   PRECONDITION(j > 0)(j <= JMax[EAIndex(NLower, 1)]);
   PRECONDITION_EQUAL(UnitCellSize, BDeque.size());

   EAIndex Index(NLower, j);

   int n = 0;
   for (auto const& B : BDeque)
      WindowLower[Index][n++] = B;

   // Invalidate elements which depend on this window element.
   for (auto I = EMatK.begin(); I != EMatK.end();)
   {
      bool Erase = true;
      for (auto const& J : I->first.second)
      {
         if (J < j)
         {
            Erase = false;
            break;
         }
      }

      if (Erase)
         I = EMatK.erase(I);
      else
         ++I;
   }

   for (auto I = FMatK.begin(); I != FMatK.end();)
   {
      bool Erase = true;
      for (auto const& J : I->first.second)
      {
         if (J >= j)
         {
            Erase = false;
            break;
         }
      }

      if (Erase)
         I = FMatK.erase(I);
      else
         ++I;
   }
}

StateComponent
EFMatrix::GetWUpper(EAIndex i, int n)
{
   // If not defined, return null. (TODO: Is this the best behavior?)
   if (WindowUpper.count(i) == 0)
      return StateComponent();

   // Get the site number of the element to convert the window position to the
   // unit cell position.
   int Site = 0;
   for (auto const& I : i)
   {
      if (!IsCorner(I))
      {
         Site = I.value();
         break;
      }
   }

   StateComponent WUpper = WindowUpper[i][numerics::divp(n-Site+1,UnitCellSize).rem];
   for (int I = numerics::divp(n-Site+1,UnitCellSize).quot; I < 0; ++I)
      WUpper.delta_shift(QShift);
   return WUpper;
}

StateComponent
EFMatrix::GetWLower(EAIndex j, int n)
{
   // If not defined, return null.
   if (WindowLower.count(j) == 0)
      return StateComponent();

   // Get the site number of the element to convert the window position to the
   // unit cell position.
   int Site = 0;
   for (auto const& I : j)
   {
      if (!IsCorner(I))
      {
         Site = I.value();
         break;
      }
   }

   StateComponent WLower = WindowLower[j][numerics::divp(n-Site+1,UnitCellSize).rem];
   for (int I = numerics::divp(n-Site+1,UnitCellSize).quot; I < 0; ++I)
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
EFMatrix::GetTLeft(CornerIndex i, CornerIndex j, bool F)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(j.size(), NLower);

   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   if (F)
      return delta_shift(TLeft[std::make_pair(i, j)], adjoint(QShift));
   else
      return TLeft[std::make_pair(i, j)];
}

MatrixOperator
EFMatrix::GetTRight(CornerIndex i, CornerIndex j, bool F)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(j.size(), NLower);

   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   if (F)
      return delta_shift(TRight[std::make_pair(i, j)], adjoint(QShift));
   else
      return TRight[std::make_pair(i, j)];
}

void
EFMatrix::CalculateTEVs(CornerIndex i, CornerIndex j)
{
   ProductMPO StringOp = ProductMPO::make_identity(ExtractLocalBasis(PsiUpper[i]));
   if (Op.is_product())
      StringOp = Op.as_product_mpo();

   // Find the eigenpair if they have an eigenvalue of magnitude 1.
   std::tie(std::ignore, TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)])
      = get_transfer_unit_eigenpair(PsiUpper[i], PsiLower[j], QShift, StringOp, Tol, UnityEpsilon, Verbose-1);

   TRight[std::make_pair(i, j)].delta_shift(QShift);

   // Special normalization for the first transfer matrix.
   if (i == CornerIndex(NUpper, 0) && j == CornerIndex(NLower, 0))
      Normalize(TRight[std::make_pair(i, j)], TLeft[std::make_pair(i, j)]);

   // Special normalization for the final transfer matrix.
   if (i == CornerIndex(NUpper, Infinity) && j == CornerIndex(NLower, Infinity))
      Normalize(TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)]);

   TCalculated[std::make_pair(i, j)] = true;
}

std::complex<double>
EFMatrix::MomentumFactor(CornerIndex i, CornerIndex j)
{
   return ExpIKUpper[i] * std::conj(ExpIKLower[j]);
}

std::complex<double>
EFMatrix::MomentumFactor(EAIndex i, EAIndex j)
{
   return this->MomentumFactor(ToCorner(i), ToCorner(j));
}

std::vector<KMatrixPolyType>
EFMatrix::GetE(EAIndex i, EAIndex j, int n)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(j.size(), NLower);

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

std::map<EFMatrix::EAIndex, StateComponent>
EFMatrix::GetWPrev(EAIndex i, int n, bool Upper)
{
   std::map<EAIndex, StateComponent> Result;

   if (!IsCorner(i))
   {
      // Get the index for the previous element
      EAIndex Index = i;
      for (auto& I : Index)
         if (!IsCorner(I))
            I = I.value() - 1;

      Result[Index] = this->GetW(i, n, Upper);
   }
   else
   {
      // For a corner element, we may have multiple windows, so we get the
      // index for the diagonal window, and construct the rest recursively.
      EAIndex DiagIndex = i;
      for (auto& I : DiagIndex)
         if (I.is_inf())
            I = 1;

      this->GetWPrevCorner(Result, DiagIndex, n, Upper);
   }

   return Result;
}

void
EFMatrix::GetWPrevCorner(std::map<EAIndex, StateComponent>& Result, EAIndex i, int n, bool Upper)
{
   // Call this function replacing each 1 in the input with infinity.
   for (auto& I : i)
   {
      if (I == 1)
      {
         I = Infinity;
         this->GetWPrevCorner(Result, i, n, Upper);
         I = 1;
      }
   }

   // Get the element for this window if it exists.
   if (Upper ? IMax.count(i) : JMax.count(i))
   {
      int Max = Upper ? IMax[i] : JMax[i];

      EAIndex Index = i;
      for (auto& I : Index)
         if (I == 1)
            I = Max;

      EAIndex IndexMinusOne = i;
      for (auto& I : IndexMinusOne)
         if (I == 1)
            I = Max-1;

      Result[IndexMinusOne] = this->GetW(Index, n, Upper);
   }
}

void
EFMatrix::CalculateE(EAIndex i, EAIndex j)
{
   auto O = Op.as_generic_mpo().begin();
   LinearWavefunction::const_iterator CUpper, CLower;

   bool CornerUpper = IsCorner(i);
   if (CornerUpper)
      CUpper = PsiUpper[ToCorner(i)].begin();

   bool CornerLower = IsCorner(j);
   if (CornerLower)
      CLower = PsiLower[ToCorner(j)].begin();

   // Cumulative sum for corner elements.
   std::vector<KMatrixPolyType> CTriK = std::vector<KMatrixPolyType>(O->Basis2().size());

   // Loop over each position in the unit cell.
   for (int n = 0; n < UnitCellSize; ++n)
   {
      EMatK[std::make_pair(i, j)][n] = std::vector<KMatrixPolyType>(O->Basis2().size());

      // Handle contributions with a window on the top and the bottom.
      for (auto const& I : this->GetWPrev(i, n, true))
         for (auto const& J : this->GetWPrev(j, n, false))
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(I.second), this->GetE(I.first, J.first, n-1), J.second);

      // Handle contributions with a window on the bottom only.
      if (CornerUpper)
         for (auto const& J : this->GetWPrev(j, n, false))
            EMatK[std::make_pair(i, j)][n] += contract_from_left(*O, herm(*CUpper), this->GetE(i, J.first, n-1), J.second);

      // Handle contributions with a window on the top only.
      if (CornerLower)
         for (auto const& I : this->GetWPrev(i, n, true))
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
      this->SolveE(ToCorner(i), ToCorner(j), delta_shift(CTriK, QShift));
}

void
EFMatrix::SolveE(CornerIndex i, CornerIndex j, std::vector<KMatrixPolyType> CTriK)
{
   std::pair<EAIndex, EAIndex> Index(EAIndex(i.begin(), i.end()), EAIndex(j.begin(), j.end()));

   EMatK[Index][UnitCellSize-1].clear();

   bool FirstElement = i == CornerIndex(NUpper, 0) && j == CornerIndex(NLower, 0);
   bool FinalElement = i == CornerIndex(NUpper, Infinity) && j == CornerIndex(NLower, Infinity);

   // Initialize the first element of the first E matrix.
   if (FirstElement)
   {
      EMatK[Index][UnitCellSize-1].push_back(KMatrixPolyType());
      EMatK[Index][UnitCellSize-1][0][1.0] = MatrixPolyType(this->GetTLeft(i, j));
   }

   // Run the linear solver.
   SolveMPO_EA_Left(EMatK[Index][UnitCellSize-1], CTriK, PsiUpper[i], PsiLower[j], QShift,
                    Op, this->GetTLeft(i, j), this->GetTRight(i, j), this->MomentumFactor(i, j),
                    Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                    !FirstElement && EAOptimization, Verbose-1);

   EMatK[Index][UnitCellSize-1] = delta_shift(EMatK[Index][UnitCellSize-1], adjoint(QShift));

   // Shift polynomial variable from n-1 to n.
   // TODO: Can this be factored out?
   ShiftVariable(EMatK[Index][UnitCellSize-1], -1);

   // Calculate the elements in the rest of the unit cell.
   if (!FinalElement || NeedFinalMatrix)
   {
      auto CUpper = PsiUpper[i].begin();
      auto CLower = PsiLower[j].begin();
      auto O = Op.as_generic_mpo().begin();

      for (int n = 0; n < UnitCellSize-1; ++n)
      {
         EMatK[Index][n] += contract_from_left(*O, herm(*CUpper), this->GetE(Index.first, Index.second, n-1), *CLower);
         ++CUpper, ++CLower, ++O;
      }
   }
}

std::vector<KMatrixPolyType>
EFMatrix::GetF(EAIndex i, EAIndex j, int n)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(j.size(), NLower);

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

std::map<EFMatrix::EAIndex, StateComponent>
EFMatrix::GetWNext(EAIndex i, int n, bool Upper)
{
   std::map<EAIndex, StateComponent> Result;

   if (!IsCorner(i))
   {
      // Get the index for the first element in this window.
      EAIndex FirstIndex = i;
      int Site = 0;
      for (auto& I : FirstIndex)
      {
         if (!IsCorner(I))
         {
            Site = I.value();
            I = 1;
         }
      }

      int Max = Upper ? IMax[FirstIndex] : JMax[FirstIndex];

      // Get the index for the next element and the next window.
      EAIndex Index = i;
      for (auto& I : Index)
         if (!IsCorner(I))
            I = Site == Max-1 ? ExtendedInt(Infinity) : ExtendedInt(Site+1);

      EAIndex WindowIndex = i;
      for (auto& I : WindowIndex)
         if (!IsCorner(I))
            I = Site+1;

      Result[Index] = this->GetW(WindowIndex, n, Upper);
   }
   else
   {
      // For a corner element, we may have multiple windows, so we get the
      // index for the diagonal window, and construct the rest recursively.
      EAIndex DiagIndex = i;
      for (auto& I : DiagIndex)
         if (I == 0)
            I = 1;

      this->GetWNextCorner(Result, DiagIndex, n, Upper);
   }

   return Result;
}

void
EFMatrix::GetWNextCorner(std::map<EAIndex, StateComponent>& Result, EAIndex i, int n, bool Upper)
{
   // Call this function replacing each 1 in the input with 0.
   for (auto& I : i)
   {
      if (I == 1)
      {
         I = 0;
         this->GetWNextCorner(Result, i, n, Upper);
         I = 1;
      }
   }

   // Get the element for this window if it exists.
   if (Upper ? IMax.count(i) : JMax.count(i))
   {
      int Max = Upper ? IMax[i] : JMax[i];

      // The index for the next element is infinity if the window size is one.
      EAIndex Index = i;
      if (Max == 1)
         for (auto& I : Index)
            if (I == 1)
               I = Infinity;

      Result[Index] = this->GetW(i, n, Upper);
   }
}

void
EFMatrix::CalculateF(EAIndex i, EAIndex j)
{
   auto O = Op.as_generic_mpo().end();
   LinearWavefunction::const_iterator CUpper, CLower;

   bool CornerUpper = IsCorner(i);
   if (CornerUpper)
      CUpper = PsiUpper[ToCorner(i)].end();

   bool CornerLower = IsCorner(j);
   if (CornerLower)
      CLower = PsiLower[ToCorner(j)].end();

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
      for (auto const& I : this->GetWNext(i, n, true))
         for (auto const& J : this->GetWNext(j, n, false))
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), I.second, this->GetF(I.first, J.first, n+1), herm(J.second));

      // Handle contributions with a window on the bottom only.
      if (CornerUpper)
         for (auto const& J : this->GetWNext(j, n, false))
            FMatK[std::make_pair(i, j)][n] += contract_from_right(herm(*O), *CUpper, this->GetF(i, J.first, n+1), herm(J.second));

      // Handle contributions with a window on the top only.
      if (CornerLower)
         for (auto const& I : this->GetWNext(i, n, true))
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
      this->SolveF(ToCorner(i), ToCorner(j), delta_shift(CTriK, adjoint(QShift)));
}

void
EFMatrix::SolveF(CornerIndex i, CornerIndex j, std::vector<KMatrixPolyType> CTriK)
{
   std::pair<EAIndex, EAIndex> Index(EAIndex(i.begin(), i.end()), EAIndex(j.begin(), j.end()));

   FMatK[Index][0].clear();

   bool FirstElement = i == CornerIndex(NUpper, Infinity) && j == CornerIndex(NLower, Infinity);
   bool FinalElement = i == CornerIndex(NUpper, 0) && j == CornerIndex(NLower, 0);

   // Initialize the first element of the first F matrix.
   if (FirstElement)
   {
      FMatK[Index][0].push_back(KMatrixPolyType());
      FMatK[Index][0][0][1.0] = MatrixPolyType(this->GetTRight(i, j, true));
   }

   // Run the linear solver.
   SolveMPO_EA_Right(FMatK[Index][0], CTriK, PsiUpper[i], PsiLower[j], QShift,
                     Op, this->GetTLeft(i, j, true), this->GetTRight(i, j, true), std::conj(this->MomentumFactor(i, j)),
                     Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                     !FirstElement && EAOptimization, Verbose-1);

   FMatK[Index][0] = delta_shift(FMatK[Index][0], QShift);

   // Subtract the contribution due to the energy density: this is to ensure
   // that the expectation value obtained by taking product of E and F matrices
   // corresponds to the value we would obtain by solving the full E or F
   // matrix.
   if (FirstElement && SubtractEnergy)
   {
      // The first contribution can be extracted from the expectation value of the F matrix.
      std::complex<double> RightEnergy = inner_prod(this->GetTLeft(i, j), FMatK[Index][0].front()[1.0].coefficient(1));

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
         BondEnergy += inner_prod(EMatKRight[I][1.0].coefficient(0), FMatK[Index][0][I][1.0].coefficient(0));

      if (Verbose > 0)
         std::cerr << "EFMatrix: subtracting contribution from bond energy = " << BondEnergy << std::endl;

      // Subtract the right and bond energies from the final element from the first F matrix.
      FMatK[Index][0].front() -= (RightEnergy + BondEnergy) * FMatK[Index][0].back();
   }

   // Shift polynomial variable from n-1 to n.
   // TODO: Can this be factored out?
   ShiftVariable(FMatK[Index][0], -1);

   // Calculate the elements in the rest of the unit cell.
   if (!FinalElement || NeedFinalMatrix)
   {
      auto CUpper = PsiUpper[i].end();
      auto CLower = PsiLower[j].end();
      auto O = Op.as_generic_mpo().end();

      for (int n = UnitCellSize-1; n > 0; --n)
      {
         --CUpper, --CLower, --O;
         FMatK[Index][n] += contract_from_right(herm(*O), *CUpper, this->GetF(Index.first, Index.second, n+1), herm(*CLower));
      }
   }
}

// TODO: Could these be implemented in a better way?
StateComponent
EFMatrix::GetESC(EAIndex i, EAIndex j, int n)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(j.size(), NLower);

   std::vector<KMatrixPolyType> E = this->GetE(i, j, n);
   MatrixOperator Tmp = E[0][1.0].coefficient(0);
   OperatorComponent O = Op.as_generic_mpo()[numerics::divp(n, UnitCellSize).rem];
   StateComponent Result(O.Basis2(), Tmp.Basis1(), Tmp.Basis2());

   for (int I = 0; I < Result.size(); ++I)
      Result[I] = E[I][1.0].coefficient(0);

   return Result;
}

StateComponent
EFMatrix::GetFSC(EAIndex i, EAIndex j, int n)
{
   PRECONDITION_EQUAL(i.size(), NUpper);
   PRECONDITION_EQUAL(j.size(), NLower);

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

   int iMax = IMax[EAIndex(NUpper, 1)];
   int jMax = JMax[EAIndex(NLower, 1)];

   auto CLeft = PsiLower[CornerIndex(NLower, 0)].begin();
   auto CRight = PsiLower[CornerIndex(NLower, Infinity)].begin();
   auto O = Op.as_generic_mpo().begin();

   ExtendedInt iNext;
   if (i == iMax-1)
      iNext = Infinity;
   else
      iNext = i+1;

   // Loop over each position in the unit cell.
   for (int n = 0; n < UnitCellSize; ++n)
   {
      // Add each contribution due to each component on the bottom.
      Result.push_back(operator_prod_inner(*O, this->GetESC({i}, {0}, n-1), *CLeft, herm(this->GetFSC({iNext}, {0}, n+1))));

      for (int j = 0; j < jMax; ++j)
      {
         ExtendedInt jNext;
         if (j == jMax-1)
            jNext = Infinity;
         else
            jNext = j+1;

         Result.back() += operator_prod_inner(*O, this->GetESC({i}, {j}, n-1), this->GetWLower({j+1}, n), herm(this->GetFSC({iNext}, {jNext}, n+1)));
      }

      Result.back() += operator_prod_inner(*O, this->GetESC({i}, {Infinity}, n-1), *CRight, herm(this->GetFSC({iNext}, {Infinity}, n+1)));

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
