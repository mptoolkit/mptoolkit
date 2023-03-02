// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/triangular_mpo_solver_ea.cpp
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "transfer.h"
#include "triangular_mpo_solver.h"
#include "triangular_mpo_solver_helpers.h"

void
SolveMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK, std::vector<KMatrixPolyType> const& CTriK,
                 LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                 QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                 MatrixOperator const& TLeft, MatrixOperator const& TRight,
                 std::complex<double> ExpIK, int Degree, double Tol, double UnityEpsilon,
                 bool NeedFinalMatrix, bool EAOptimization, int Verbose)
{
   int Dim = Op.Basis1().size();
   int StartCol = EMatK.size();
   EMatK.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_EA_Left: dimension = " << Dim << std::endl;

   // Solve recursively from the first empty column.
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
         std::cerr << "Solving column " << (Col) << " of [0:" << (Dim-1) << "]" << std::endl;

      // Generate the next C matrices.
      KMatrixPolyType C = ExpIK * inject_left_mask(EMatK, Psi1, QShift, Op.data(), Psi2, mask_column(Op, Col))[Col];

      // Add terms upper-triangular in the EA indices.
      if (!CTriK.empty())
         C += CTriK[Col];

      // Now do the classification, based on the properties of the diagonal operator.
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;

         EMatK[Col] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

         // Components parallel to the transfer matrix left eigenvector.
         KComplexPolyType EParallel;

         // Multiplication factor and left and right transfer matrix eigenvectors.
         std::complex<double> Factor = Classification.factor();
         MatrixOperator TransferEVLeft = TLeft;
         MatrixOperator TransferEVRight = TRight;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            // Find the eigenvectors of the transfer matrix: if the operator is
            // proportional to the identity in the scalar sector, then the
            // eigenvectors are those of the usual transfer matrix (TLeft and
            // TRight), otherwise, we must solve for the generalised transfer
            // matrix eigenvalues.
            if (!is_scalar(Diag.Basis2()[0]) || !Classification.is_complex_identity())
            {
               if (Verbose > 0)
                  std::cerr << "Solving unitary diagonal component" << std::endl;

               // Find the largest eigenvalue.
               std::complex<double> EValue;
               std::tie(EValue, TransferEVLeft, TransferEVRight)
                  = get_transfer_unit_eigenpair(Psi1, Psi2, QShift, ProductMPO(Diag), Tol, UnityEpsilon, Verbose);
               TransferEVRight = delta_shift(TransferEVRight, QShift);

               EValue = std::conj(EValue); // left eigenvalue, so conjugate (see comment at operator_actions.h)

               if (Verbose > 0)
                  std::cerr << "Eigenvalue of unitary operator is " << EValue << std::endl;

               Factor = EValue;
            }
            else if (Verbose > 0 && Classification.is_identity())
               std::cerr << "Diagonal component is the identity" << std::endl;
            else if (Verbose > 0 && Classification.is_complex_identity())
               std::cerr << "Diagonal component is proportional to the identity" << std::endl;

            // Multiply by the momentum factor.
            Factor *= ExpIK;

            // We only need to solve for the parallel components if we have an
            // eigenvalue of magnitude 1, which we can tell by whether the
            // eigenvector is null.
            if (!TransferEVLeft.is_null())
            {
               HasEigenvalue1 = true;

               if (Verbose > 0)
                  std::cerr << "Decomposing parts parallel to the unit matrix" << std::endl;

               EParallel = DecomposeParallelPartsWithMomentum(C, Factor, TransferEVLeft, TransferEVRight, UnityEpsilon, Degree);
            }
            else if (Verbose > 0)
               std::cerr << "Diagonal component has spectral radius < 1" << std::endl;
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1.

         KMatrixPolyType E;
         // If we are on the last column and we don't need the matrix elements, then we can
         // skip this operation.
         if (Col < Dim-1 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix" << std::endl;

            if (EAOptimization && Col == 0)
               E[1.0] = 0.0 * C[1.0]; // This will be zero if we are in the left gauge.
            else if (EAOptimization && Col == Dim-1)
            {
               // This should be zero anyway, so we do not want to run the linear solver for this component.
               C[1.0].erase(1);
               // For the EA algorithm, we only need the zero momentum components for the final column.
               E[1.0] = DecomposePerpendicularPartsLeft(C[1.0], 1.0, ExpIK*Diag, TransferEVLeft, TransferEVRight,
                                                        Psi1, Psi2, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            }
            else
               E = DecomposePerpendicularPartsLeft(C, ExpIK*Diag, TransferEVLeft, TransferEVRight,
                                                   Psi1, Psi2, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last column" << std::endl;

         // Reinsert the components parallel to the unit matrix (if any).
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               // Conj here because this comes from an overlap(x, TransferEVRight).
               E[I->first][J->first] += std::conj(J->second) * TransferEVLeft;
            }
         }

         // Finally, set the E matrix element at this column.
         EMatK[Col] = E;
      }
   }
}

void
SolveMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK, std::vector<KMatrixPolyType> const& CTriK,
                  LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                  QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                  MatrixOperator const& TLeft, MatrixOperator const& TRight,
                  std::complex<double> ExpIK, int Degree, double Tol, double UnityEpsilon,
                  bool NeedFinalMatrix, bool EAOptimization, int Verbose)
{
   int Dim = Op.Basis1().size();
   int StartRow = Dim-1-FMatK.size();
   CHECK(StartRow >= -1);
   std::vector<KMatrixPolyType> Tmp = FMatK;
   FMatK = std::vector<KMatrixPolyType>(StartRow+1);
   FMatK.insert(FMatK.end(), Tmp.begin(), Tmp.end());

   if (Verbose > 0)
      std::cerr << "SolveMPO_EA_Right: dimension = " << Dim << std::endl;

   // Solve recursively from the last empty row.
   for (int Row = StartRow; Row >= 0; --Row)
   {
      if (Verbose > 0)
         std::cerr << "Solving row " << (Row) << " of [0:" << (Dim-1) << "]" << std::endl;

      // Generate the next C matrices.
      KMatrixPolyType C = ExpIK * inject_right_mask(FMatK, Psi1, QShift, Op.data(), Psi2, mask_row(Op, Row))[Row];

      // Add terms upper-triangular in the EA indices.
      if (!CTriK.empty())
         C += CTriK[Row];

      // Now do the classification, based on the properties of the diagonal operator.
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;

         FMatK[Row] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at row " << Row << std::endl;

         // Components parallel to the transfer matrix right eigenvector.
         KComplexPolyType FParallel;

         // Multiplication factor and left and right transfer matrix eigenvectors.
         std::complex<double> Factor = Classification.factor();
         MatrixOperator TransferEVLeft = TLeft;
         MatrixOperator TransferEVRight = TRight;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            // Find the eigenvectors of the transfer matrix: if the operator is
            // proportional to the identity in the scalar sector, then the
            // eigenvectors are those of the usual transfer matrix (TLeft and
            // TRight), otherwise, we must solve for the generalised transfer
            // matrix eigenvalues.
            if (!is_scalar(Diag.Basis1()[0]) || !Classification.is_complex_identity())
            {
               if (Verbose > 0)
                  std::cerr << "Solving unitary diagonal component" << std::endl;

               // Find the largest eigenvalue.
               std::complex<double> EValue;
               std::tie(EValue, TransferEVLeft, TransferEVRight)
                  = get_transfer_unit_eigenpair(Psi1, Psi2, QShift, ProductMPO(Diag), Tol, UnityEpsilon, Verbose);
               TransferEVRight = delta_shift(TransferEVRight, QShift);

               //TransferEVLeft *= 1.0 / norm_frob(TransferEVLeft);
               //TransferEVRight *= 1.0 / inner_prod(TransferEVLeft, TransferEVRight);

               if (Verbose > 0)
                  std::cerr << "Eigenvalue of unitary operator is " << EValue << std::endl;

               Factor = EValue;
            }
            else if (Verbose > 0 && Classification.is_identity())
               std::cerr << "Diagonal component is the identity" << std::endl;
            else if (Verbose > 0 && Classification.is_complex_identity())
               std::cerr << "Diagonal component is proportional to the identity" << std::endl;

            // Multiply by the momentum factor.
            Factor *= ExpIK;

            // We only need to solve for the parallel components if we have an
            // eigenvalue of magnitude 1, which we can tell by whether the
            // eigenvector is null.
            if (!TransferEVLeft.is_null())
            {
               HasEigenvalue1 = true;

               if (Verbose > 0)
                  std::cerr << "Decomposing parts parallel to the unit matrix" << std::endl;

               FParallel = DecomposeParallelPartsWithMomentum(C, Factor, TransferEVRight, TransferEVLeft, UnityEpsilon, Degree);
            }
            else if (Verbose > 0)
               std::cerr << "Diagonal component has spectral radius < 1" << std::endl;
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1.

         KMatrixPolyType F;
         // If we are on the last row and we don't need the matrix elements, then we can
         // skip this operation.
         if (Row > 0 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix" << std::endl;

            if (EAOptimization && Row == 0)
            {
               // This should be zero anyway, so we do not want to run the linear solver for this component.
               C[1.0].erase(1);
               // For the EA algorithm, we only need the zero momentum components for the first row.
               F[1.0] = DecomposePerpendicularPartsRight(C[1.0], 1.0, std::conj(ExpIK)*Diag, TransferEVRight, TransferEVLeft,
                                                         Psi1, Psi2, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            }
            else
               F = DecomposePerpendicularPartsRight(C, std::conj(ExpIK)*Diag, TransferEVRight, TransferEVLeft,
                                                    Psi1, Psi2, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last row" << std::endl;

         // Reinsert the components parallel to the unit matrix (if any).
         for (KComplexPolyType::const_iterator I = FParallel.begin(); I != FParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               // Conj here because this comes from an overlap(x, TransferEVLeft).
               F[I->first][J->first] += std::conj(J->second) * TransferEVRight;
            }
         }

         // Finally, set the F matrix element at this row.
         FMatK[Row] = F;
      }
   }
}

void
SolveMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK, std::vector<KMatrixPolyType> const& CTriK,
                 LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                 QuantumNumber const& QShift, ProductMPO const& Op,
                 MatrixOperator const& TLeft, MatrixOperator const& TRight,
                 std::complex<double> ExpIK, int Degree, double Tol, double UnityEpsilon,
                 int Verbose)
{
   PRECONDITION(Op.is_string());

   // If EMatK is nonempty, we must already have the result, so we can return early.
   if (!EMatK.empty())
      return;

   if (Verbose > 0)
      std::cerr << "SolveMPO_EA_Left: string operator" << std::endl;

   KMatrixPolyType C = CTriK[0];

   // Components parallel to the transfer matrix left eigenvector.
   if (Verbose > 0)
      std::cerr << "Decomposing parts parallel to the unit matrix" << std::endl;

   KComplexPolyType EParallel = DecomposeParallelPartsWithMomentum(C, ExpIK, TLeft, TRight, UnityEpsilon, Degree);

   // Now the remaining components, which is anything that is not proportional
   // to an eigenvector of magnitude 1.
   if (Verbose > 0)
      std::cerr << "Decomposing parts perpendicular to the unit matrix" << std::endl;

   KMatrixPolyType E = DecomposePerpendicularPartsLeft(C, ExpIK*BasicFiniteMPO(Op), TLeft, TRight,
                                                       Psi1, Psi2, QShift, 1.0, true, Tol, Verbose);

   // Reinsert the components parallel to the unit matrix (if any).
   for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
   {
      for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
      {
         // Conj here because this comes from an overlap(x, TRight).
         E[I->first][J->first] += std::conj(J->second) * TLeft;
      }
   }

   EMatK.push_back(E);
}

void
SolveMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK, std::vector<KMatrixPolyType> const& CTriK,
                  LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                  QuantumNumber const& QShift, ProductMPO const& Op,
                  MatrixOperator const& TLeft, MatrixOperator const& TRight,
                  std::complex<double> ExpIK, int Degree, double Tol, double UnityEpsilon,
                  int Verbose)
{
   PRECONDITION(Op.is_string());

   // If FMatK is nonempty, we must already have the result, so we can return early.
   if (!FMatK.empty())
      return;

   if (Verbose > 0)
      std::cerr << "SolveMPO_EA_Right: string operator" << std::endl;

   KMatrixPolyType C = CTriK[0];

   // Components parallel to the transfer matrix right eigenvector.
   if (Verbose > 0)
      std::cerr << "Decomposing parts parallel to the unit matrix" << std::endl;

   KComplexPolyType FParallel = DecomposeParallelPartsWithMomentum(C, ExpIK, TRight, TLeft, UnityEpsilon, Degree);

   // Now the remaining components, which is anything that is not proportional
   // to an eigenvector of magnitude 1.
   if (Verbose > 0)
      std::cerr << "Decomposing parts perpendicular to the unit matrix" << std::endl;

   KMatrixPolyType F = DecomposePerpendicularPartsRight(C, std::conj(ExpIK)*BasicFiniteMPO(Op), TRight, TLeft,
                                                        Psi1, Psi2, QShift, 1.0, true, Tol, Verbose);

   // Reinsert the components parallel to the unit matrix (if any).
   for (KComplexPolyType::const_iterator I = FParallel.begin(); I != FParallel.end(); ++I)
   {
      for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
      {
         // Conj here because this comes from an overlap(x, TLeft).
         F[I->first][J->first] += std::conj(J->second) * TRight;
      }
   }

   FMatK.push_back(F);
}

void
SolveMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK, std::vector<KMatrixPolyType> const& CTriK,
                 LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                 QuantumNumber const& QShift, InfiniteMPO const& Op,
                 MatrixOperator const& TLeft, MatrixOperator const& TRight,
                 std::complex<double> ExpIK, int Degree, double Tol, double UnityEpsilon,
                 bool NeedFinalMatrix, bool EAOptimization, int Verbose)
{
   if (Op.is_triangular())
      SolveMPO_EA_Left(EMatK, CTriK, Psi1, Psi2, QShift, Op.as_basic_triangular_mpo(),
                       TLeft, TRight, ExpIK, Degree, Tol, UnityEpsilon,
                       NeedFinalMatrix, EAOptimization, Verbose);
   else if (Op.is_product())
      SolveMPO_EA_Left(EMatK, CTriK, Psi1, Psi2, QShift, Op.as_product_mpo(),
                       TLeft, TRight, ExpIK, Degree, Tol, UnityEpsilon, Verbose);
   else
      throw std::runtime_error("SolveMPO_EA_Left: fatal: InfiniteMPO must be a triangular or product MPO.");
}

void
SolveMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK, std::vector<KMatrixPolyType> const& CTriK,
                  LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                  QuantumNumber const& QShift, InfiniteMPO const& Op,
                  MatrixOperator const& TLeft, MatrixOperator const& TRight,
                  std::complex<double> ExpIK, int Degree, double Tol, double UnityEpsilon,
                  bool NeedFinalMatrix, bool EAOptimization, int Verbose)
{
   if (Op.is_triangular())
      SolveMPO_EA_Right(FMatK, CTriK, Psi1, Psi2, QShift, Op.as_basic_triangular_mpo(),
                        TLeft, TRight, ExpIK, Degree, Tol, UnityEpsilon,
                        NeedFinalMatrix, EAOptimization, Verbose);
   else if (Op.is_product())
      SolveMPO_EA_Right(FMatK, CTriK, Psi1, Psi2, QShift, Op.as_product_mpo(),
                        TLeft, TRight, ExpIK, Degree, Tol, UnityEpsilon, Verbose);
   else
      throw std::runtime_error("SolveMPO_EA_Right: fatal: InfiniteMPO must be a triangular or product MPO.");
}
