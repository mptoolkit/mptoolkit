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

// Get the window unit cell for a deque of single-site windows between two
// boundaries PsiLeft and PsiRight.
LinearWavefunction
ConstructPsiTri(LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                std::deque<StateComponent> const& BDeque)
{
   LinearWavefunction PsiTri;

   if (PsiLeft.size() == 1)
      PsiTri.push_back(BDeque.back());
   else
   {
      auto CL = PsiLeft.begin();
      auto CR = PsiRight.begin();
      auto B = BDeque.begin();
      SumBasis<VectorBasis> NewBasis0((*CL).Basis2(), (*B).Basis2());
      PsiTri.push_back(tensor_row_sum(*CL, *B, NewBasis0));
      ++CL, ++CR, ++B;
      for (int i = 1; i < PsiLeft.size()-1; ++i)
      {
         StateComponent Z = StateComponent((*CL).LocalBasis(), (*CR).Basis1(), (*CL).Basis2());
         SumBasis<VectorBasis> NewBasis1((*CL).Basis2(), (*B).Basis2());
         SumBasis<VectorBasis> NewBasis2((*CL).Basis1(), (*CR).Basis1());
         PsiTri.push_back(tensor_col_sum(tensor_row_sum(*CL, *B, NewBasis1), tensor_row_sum(Z, *CR, NewBasis1), NewBasis2));
         ++CL, ++CR, ++B;
      }
      SumBasis<VectorBasis> NewBasis3((*B).Basis1(), (*CR).Basis1());
      PsiTri.push_back(tensor_col_sum(*B, *CR, NewBasis3));
   }

   return PsiTri;
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

   // Since this is a pure virtual function, we must call it in the
   // constructors for the derived classes.
   //this->CheckOperator();
}

void
EMatrix::CheckOperator()
{
   if (Op.is_product())
   {
      // We only handle string operators at the moment.
      if (!Op.as_product_mpo().is_string())
         throw std::runtime_error("EMatrix: fatal: Op cannot be a non-string ProductMPO.");
      // We only handle scalar operators at the moment.
      if (!Op.as_product_mpo().is_scalar())
         throw std::runtime_error("EMatrix: fatal: Op cannot be a non-scalar ProductMPO.");
   }
   else if (Op.is_triangular())
   {
      if (!classify(Op.as_basic_triangular_mpo()(0,0), UnityEpsilon).is_identity())
         throw std::runtime_error("EMatrix: fatal: first component of Op must be the identity operator.");
   }
   else
      throw std::runtime_error("EMatrix: fatal: Op must be a BasicTriangularMPO or a ProductMPO.");
}

void
FMatrix::CheckOperator()
{
   if (Op.is_product())
   {
      // We only handle string operators at the moment.
      if (!Op.as_product_mpo().is_string())
         throw std::runtime_error("FMatrix: fatal: Op cannot be a non-string ProductMPO.");
      // We only handle scalar operators at the moment.
      if (!Op.as_product_mpo().is_scalar())
         throw std::runtime_error("FMatrix: fatal: Op cannot be a non-scalar ProductMPO.");
   }
   else if (Op.is_triangular())
   {
      int Dim = Op.as_basic_triangular_mpo().Basis1().size();
      if (!classify(Op.as_basic_triangular_mpo()(Dim-1,Dim-1), UnityEpsilon).is_identity())
         throw std::runtime_error("FMatrix: fatal: final component of Op must be the identity operator.");
   }
   else
      throw std::runtime_error("FMatrix: fatal: Op must be a BasicTriangularMPO or a ProductMPO.");
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
EMatrix::SetDiagTEVsLC(int i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TRight[std::make_pair(i, i)] = delta_shift(Lambda*Lambda, QShift);
   TCalculated[std::make_pair(i, i)] = true;
}

void
FMatrix::SetDiagTEVsLC(int i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis2());
   TRight[std::make_pair(i, i)] = Lambda*Lambda;
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
EMatrix::SetDiagTEVsRC(int i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = Lambda*Lambda;
   TRight[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis1());
   TCalculated[std::make_pair(i, i)] = true;
}

void
FMatrix::SetDiagTEVsRC(int i, RealDiagonalOperator Lambda)
{
   TLeft[std::make_pair(i, i)] = delta_shift(Lambda*Lambda, adjoint(QShift));
   TRight[std::make_pair(i, i)] = MatrixOperator::make_identity(PsiUpper[i].Basis2());
   TCalculated[std::make_pair(i, i)] = true;
}

void
EFMatrix::SetPsiTriUpper(int i, LinearWavefunction const& PsiTri, std::complex<double> ExpIK)
{
   PRECONDITION(i > 0);

   // Update the maximum index.
   if (i > IMax)
      IMax = i;

   PsiTriUpper[i] = PsiTri;

   // We must have already set i-1.
   // FIXME: Is it possible to do this in another way?
   CHECK(ExpIKUpper.count(i-1) == 1);
   ExpIKUpper[i] = ExpIKUpper[i-1] * ExpIK;

   // TODO: Invalidation.
}

void
EFMatrix::SetPsiTriLower(int j, LinearWavefunction const& PsiTri, std::complex<double> ExpIK)
{
   PRECONDITION(j > 0);

   // Update the maximum index.
   if (j > JMax)
      JMax = j;

   PsiTriLower[j] = PsiTri;

   // We must have already set j-1.
   // FIXME: Is it possible to do this in another way?
   CHECK(ExpIKLower.count(j-1) == 1);
   ExpIKLower[j] = ExpIKLower[j-1] * ExpIK;
}

void
EMatrix::SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK)
{
   PRECONDITION(i > 0);

   CHECK(PsiUpper.count(i-1) == 1)(PsiUpper.count(i) == 1);
   LinearWavefunction PsiTri = ConstructPsiTri(PsiUpper[i-1], PsiUpper[i], BDeque);

   this->EFMatrix::SetPsiTriUpper(i, PsiTri, ExpIK);
}

void
EMatrix::SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK)
{
   PRECONDITION(j > 0);

   CHECK(PsiLower.count(j-1) == 1)(PsiLower.count(j) == 1);
   LinearWavefunction PsiTri = ConstructPsiTri(PsiLower[j-1], PsiLower[j], BDeque);

   this->EFMatrix::SetPsiTriLower(j, PsiTri, ExpIK);
}

void
FMatrix::SetPsiTriUpper(int i, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK)
{
   PRECONDITION(i > 0);

   CHECK(PsiUpper.count(i-1) == 1)(PsiUpper.count(i) == 1);
   LinearWavefunction PsiTri = ConstructPsiTri(PsiUpper[i], PsiUpper[i-1], BDeque);

   this->EFMatrix::SetPsiTriUpper(i, PsiTri, ExpIK);
}

void
FMatrix::SetPsiTriLower(int j, std::deque<StateComponent> const& BDeque, std::complex<double> ExpIK)
{
   PRECONDITION(j > 0);

   CHECK(PsiLower.count(j-1) == 1)(PsiLower.count(j) == 1);
   LinearWavefunction PsiTri = ConstructPsiTri(PsiLower[j], PsiLower[j-1], BDeque);

   this->EFMatrix::SetPsiTriLower(j, PsiTri, ExpIK);
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
   EFMatK.clear();

   Op = Op_;
   Degree = Degree_;

   this->CheckOperator();
}

MatrixOperator
EFMatrix::GetTLeft(int i, int j)
{
   // Return null if out of bounds.
   if (i < 0 || j < 0)
      return MatrixOperator();

   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   return TLeft[std::make_pair(i, j)];
}

MatrixOperator
EFMatrix::GetTRight(int i, int j)
{
   // Return null if out of bounds.
   if (i < 0 || j < 0)
      return MatrixOperator();

   // Calculate the element if we haven't already.
   if (!TCalculated[std::make_pair(i, j)])
      this->CalculateTEVs(i, j);

   return TRight[std::make_pair(i, j)];
}

void
EMatrix::CalculateTEVs(int i, int j)
{
   ProductMPO StringOp = ProductMPO::make_identity(ExtractLocalBasis(PsiUpper[i]));
   if (Op.is_product())
      StringOp = Op.as_product_mpo();

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

void
FMatrix::CalculateTEVs(int i, int j)
{
   ProductMPO StringOp = ProductMPO::make_identity(ExtractLocalBasis(PsiUpper[i]));
   if (Op.is_product())
      StringOp = Op.as_product_mpo();

   std::tie(std::ignore, TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)])
      = get_transfer_unit_eigenpair(PsiUpper[i], PsiLower[j], QShift, StringOp, Tol, UnityEpsilon, Verbose);

   TLeft[std::make_pair(i, j)].delta_shift(adjoint(QShift));

   // Special normalization for the first transfer matrix.
   if (i == 0 && j == 0)
      Normalize(TLeft[std::make_pair(i, j)], TRight[std::make_pair(i, j)]);

   // Special normalization for the final transfer matrix.
   if (i == IMax && j == JMax)
      Normalize(TRight[std::make_pair(i, j)], TLeft[std::make_pair(i, j)]);

   TCalculated[std::make_pair(i, j)] = true;
}

std::vector<KMatrixPolyType>
EFMatrix::GetElement(int i, int j)
{
   // Return null if out of bounds.
   if (i < 0 || j < 0)
      return std::vector<KMatrixPolyType>();

   // Calculate the element if we haven't already.
   if (EFMatK[std::make_pair(i, j)].empty())
      this->CalculateElement(i, j);

   return EFMatK[std::make_pair(i, j)];
}

void
EMatrix::CalculateElement(int i, int j)
{
   // Initialize the first element of the first E matrix.
   if (i == 0 && j == 0)
   {
      EFMatK[std::make_pair(i, j)].push_back(KMatrixPolyType());
      EFMatK[std::make_pair(i, j)][0][1.0] = MatrixPolyType(this->GetIdent(i, j));
   }

   std::vector<KMatrixPolyType> CTriK
      = CalculateCTriK_Left(this->GetElement(i, j-1), this->GetElement(i-1, j), this->GetElement(i-1, j-1),
                            PsiUpper[i], PsiLower[j], PsiTriUpper[i], PsiTriLower[j],
                            QShift, Op.as_generic_mpo(), ExpIKUpper[i], ExpIKLower[j]);

   bool FinalElement = i == IMax && j == JMax;

   SolveMPO_EA_Left(EFMatK[std::make_pair(i, j)], CTriK, PsiUpper[i], PsiLower[j],
                    QShift, Op, this->GetTLeft(i, j), this->GetTRight(i, j),
                    ExpIKUpper[i] * std::conj(ExpIKLower[j]),
                    Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                    FinalElement && EAOptimization, Verbose);
}

void
FMatrix::CalculateElement(int i, int j)
{
   // Initialize the first element of the first F matrix.
   if (i == 0 && j == 0)
   {
      EFMatK[std::make_pair(i, j)].push_back(KMatrixPolyType());
      EFMatK[std::make_pair(i, j)][0][1.0] = MatrixPolyType(this->GetIdent(i, j));
   }

   std::vector<KMatrixPolyType> CTriK
      = CalculateCTriK_Right(this->GetElement(i, j-1), this->GetElement(i-1, j), this->GetElement(i-1, j-1),
                             PsiUpper[i], PsiLower[j], PsiTriUpper[i], PsiTriLower[j],
                             QShift, Op.as_generic_mpo(), ExpIKUpper[i], ExpIKLower[j]);

   bool FinalElement = i == IMax && j == JMax;

   SolveMPO_EA_Right(EFMatK[std::make_pair(i, j)], CTriK, PsiUpper[i], PsiLower[j],
                     QShift, Op, this->GetTLeft(i, j), this->GetTRight(i, j),
                     ExpIKUpper[i] * std::conj(ExpIKLower[j]),
                     Degree, Tol, UnityEpsilon, !FinalElement || NeedFinalMatrix,
                     FinalElement && EAOptimization, Verbose);
}
