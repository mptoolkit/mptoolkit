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
SolveMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK1, std::vector<KMatrixPolyType> const& EMatK0,
                 std::vector<KMatrixPolyType> const& EMatKTop, std::vector<KMatrixPolyType> const& EMatKBot,
                 LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                 LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                 BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                 MatrixOperator const& TRight, std::complex<double> ExpIK,
                 bool NeedFinalMatrix, int Degree, double Tol,
                 double UnityEpsilon, std::string Mode, int Verbose)
{
   CHECK(Mode == "initial" || Mode == "top" || Mode == "top-ea" || Mode == "bottom" || Mode == "final");

   if (Mode == "initial" || Mode == "final")
      CHECK(!TLeft.is_null() && !TRight.is_null());

   int Dim = Op.Basis1().size();
   EMatK1.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_EA_Left: mode = " << Mode << ", dimension = " << Dim << std::endl;

   int StartCol = 0;

   if (EMatK1[0][1.0].empty() && Mode == "initial")
   {
      // Make sure the top-left element is the identity.
      if (!classify(Op(StartCol, StartCol), UnityEpsilon).is_identity())
      {
         //std::cerr << "SolveMPO_EA_Left: fatal: MPO(0,0) must be the identity operator for mode initial." << std::endl;
         //PANIC("Fatal");
      }

      // Initialize the first E-matrix.
      EMatK1[StartCol][1.0] = MatrixPolyType(TLeft);

      ++StartCol;
   }

   // Solve recursively from the first column.
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
         std::cerr << "Solving column " << (Col) << " of [0:" << (Dim-1) << "]" << std::endl;

      // Generate the next C matrices: this depends on the E-matrix being calculated.
      KMatrixPolyType C;

      if (Mode == "initial")
      {
         C = inject_left_mask(EMatK1, PsiLeft, QShift, Op.data(), PsiLeft, mask_column(Op, Col))[Col];
      }
      else if (Mode == "top" || Mode == "top-ea")
      {
         C = inject_left_mask(EMatK1, PsiRight, QShift, Op.data(), PsiLeft, mask_column(Op, Col))[Col];
         for (auto I = C.begin(); I != C.end(); ++I)
            I->second *= ExpIK;
         KMatrixPolyType CInit = inject_left_mask(EMatK0, PsiTri, QShift, Op.data(), PsiLeft, mask_column(Op, Col))[Col];
         for (auto I = CInit.begin(); I != CInit.end(); ++I)
            C[I->first] += I->second;
      }
      else if (Mode == "bottom")
      {
         C = inject_left_mask(EMatK1, PsiLeft, QShift, Op.data(), PsiRight, mask_column(Op, Col))[Col];
         for (auto I = C.begin(); I != C.end(); ++I)
            I->second *= std::conj(ExpIK);
         KMatrixPolyType CInit = inject_left_mask(EMatK0, PsiLeft, QShift, Op.data(), PsiTri, mask_column(Op, Col))[Col];
         for (auto I = CInit.begin(); I != CInit.end(); ++I)
            C[I->first] += I->second;
      }
      else if (Mode == "final")
      {
         C = inject_left_mask(EMatK1, PsiRight, QShift, Op.data(), PsiRight, mask_column(Op, Col))[Col];
         KMatrixPolyType CTop = inject_left_mask(EMatKTop, PsiRight, QShift, Op.data(), PsiTri, mask_column(Op, Col))[Col];
         for (auto I = CTop.begin(); I != CTop.end(); ++I)
            C[I->first] += ExpIK * I->second;
         KMatrixPolyType CBot = inject_left_mask(EMatKBot, PsiTri, QShift, Op.data(), PsiRight, mask_column(Op, Col))[Col];
         for (auto I = CBot.begin(); I != CBot.end(); ++I)
            C[I->first] += std::conj(ExpIK) * I->second;
         KMatrixPolyType CInit = inject_left_mask(EMatK0, PsiTri, QShift, Op.data(), PsiTri, mask_column(Op, Col))[Col];
         for (auto I = CInit.begin(); I != CInit.end(); ++I)
            C[I->first] += I->second;
      }

      // Now do the classification, based on the properties of the diagonal operator.
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;

         EMatK1[Col] = SolveZeroDiagonal(C);
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
         if (Classification.is_unitary() || Classification.is_unclassified())
         {
            // Find the eigenvectors of the transfer matrix: if the operator is
            // proportional to the identity in the scalar sector, then the
            // eigenvectors are those of the usual transfer matrix (TLeft and
            // TRight), otherwise, we must solve for the generalised transfer
            // matrix eigenvalues.
            if (!is_scalar(Diag.Basis2()[0]) || !(Classification.is_complex_identity() || Classification.is_unclassified()))
            {
               if (Verbose > 0)
                  std::cerr << "Solving unitary diagonal component" << std::endl;

               // Find the largest eigenvalue.
               std::complex<double> EValue;
               if (Mode == "initial")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiLeft, PsiLeft, QShift, ProductMPO(Diag), Tol, Verbose);
               else if (Mode == "top" || Mode == "top-ea")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiRight, PsiLeft, QShift, ProductMPO(Diag), Tol, Verbose);
               else if (Mode == "bottom")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiLeft, PsiRight, QShift, ProductMPO(Diag), Tol, Verbose);
               else if (Mode == "final")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiRight, PsiRight, QShift, ProductMPO(Diag), Tol, Verbose);

               TransferEVRight = delta_shift(TransferEVRight, QShift);

               EValue = std::conj(EValue); // left eigenvalue, so conjugate (see comment at operator_actions.h)

               if (Verbose > 0)
                  std::cerr << "Eigenvalue of unitary operator is " << EValue << std::endl;

               // If the magnitude of the largest eigenvalue is less than one, we don't need the eigenvectors.
               if (std::abs(norm_frob(Factor) - 1.0) < UnityEpsilon)
               {
                  TransferEVLeft = MatrixOperator();
                  TransferEVRight = MatrixOperator();
               }
               else
                  Factor = EValue;
            }
            else if (Verbose > 0 && Classification.is_identity())
               std::cerr << "Diagonal component is the identity" << std::endl;
            else if (Verbose > 0 && Classification.is_complex_identity())
               std::cerr << "Diagonal component is proportional to the identity" << std::endl;
            else if (Classification.is_unclassified())
            {
               if (Verbose > 0)
                  std::cerr << "Diagonal component is unclassified" << std::endl;
               Factor = 1.0;
            }

            // Multiply by the EA momentum factors for mixed transfer matrices.
            if (Mode == "top" || Mode == "top-ea")
               Factor *= ExpIK;
            else if (Mode == "bottom")
               Factor *= std::conj(ExpIK);

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

            if (Mode == "initial")
               E = DecomposePerpendicularPartsLeft(C, Diag, TransferEVLeft, TransferEVRight,
                                                   PsiLeft, PsiLeft, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "top-ea" && Col == 0 && Classification.is_identity())
               E[1.0] = 0.0 * C[1.0]; // This will be zero if we are in the left gauge.
            else if (Mode == "top-ea" && Col == Dim-1)
               // For the EA algorithm, we only need the zero momentum components for the final column.
               E[1.0] = DecomposePerpendicularPartsLeft(C[1.0], 1.0, ExpIK*Diag, TransferEVLeft, TransferEVRight,
                                                        PsiRight, PsiLeft, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "top" || Mode == "top-ea")
               E = DecomposePerpendicularPartsLeft(C, ExpIK*Diag, TransferEVLeft, TransferEVRight,
                                                   PsiRight, PsiLeft, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "bottom")
               E = DecomposePerpendicularPartsLeft(C, std::conj(ExpIK)*Diag, TransferEVLeft, TransferEVRight,
                                                   PsiLeft, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "final")
               E = DecomposePerpendicularPartsLeft(C, Diag, TransferEVLeft, TransferEVRight,
                                                   PsiRight, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
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
         EMatK1[Col] = E;
      }
   }
}

void
SolveMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK1, std::vector<KMatrixPolyType> const& FMatK0,
                  std::vector<KMatrixPolyType> const& FMatKTop, std::vector<KMatrixPolyType> const& FMatKBot,
                  LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                  LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                  BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                  MatrixOperator const& TRight, std::complex<double> ExpIK,
                  bool NeedFinalMatrix, int Degree, double Tol,
                  double UnityEpsilon, std::string Mode, int Verbose)
{
   CHECK(Mode == "initial" || Mode == "top" || Mode == "top-ea" || Mode == "bottom" || Mode == "final");

   if (Mode == "initial" || Mode == "final")
      CHECK(!TLeft.is_null() && !TRight.is_null());

   int Dim = Op.Basis1().size();
   FMatK1.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_EA_Right: mode = " << Mode << ", dimension = " << Dim << std::endl;

   int StartRow = Dim-1;

   if (FMatK1[0][1.0].empty() && Mode == "initial")
   {
      // Make sure the bottom-right element is the identity.
      if (!classify(Op(StartRow, StartRow), UnityEpsilon).is_identity())
      {
         //std::cerr << "SolveMPO_EA_Right: fatal: MPO(last,last) must be the identity operator for mode initial." << std::endl;
         //PANIC("Fatal");
      }

      // Initialize the first F-matrix.
      FMatK1[StartRow][1.0] = MatrixPolyType(TRight);

      --StartRow;
   }

   // Solve recursively from the final row.
   for (int Row = StartRow; Row >= 0; --Row)
   {
      if (Verbose > 0)
         std::cerr << "Solving row " << (Row) << " of [0:" << (Dim-1) << "]" << std::endl;

      // Generate the next C matrices: this depends on the F-matrix being calculated.
      KMatrixPolyType C;

      if (Mode == "initial")
      {
         C = inject_right_mask(FMatK1, PsiRight, QShift, Op.data(), PsiRight, mask_row(Op, Row))[Row];
      }
      else if (Mode == "top" || Mode == "top-ea")
      {
         C = inject_right_mask(FMatK1, PsiLeft, QShift, Op.data(), PsiRight, mask_row(Op, Row))[Row];
         for (auto I = C.begin(); I != C.end(); ++I)
            I->second *= ExpIK;
         KMatrixPolyType CInit = inject_right_mask(FMatK0, PsiTri, QShift, Op.data(), PsiRight, mask_row(Op, Row))[Row];
         for (auto I = CInit.begin(); I != CInit.end(); ++I)
            C[I->first] += I->second;
      }
      else if (Mode == "bottom")
      {
         C = inject_right_mask(FMatK1, PsiRight, QShift, Op.data(), PsiLeft, mask_row(Op, Row))[Row];
         for (auto I = C.begin(); I != C.end(); ++I)
            I->second *= std::conj(ExpIK);
         KMatrixPolyType CInit = inject_right_mask(FMatK0, PsiRight, QShift, Op.data(), PsiTri, mask_row(Op, Row))[Row];
         for (auto I = CInit.begin(); I != CInit.end(); ++I)
            C[I->first] += I->second;
      }
      else if (Mode == "final")
      {
         C = inject_right_mask(FMatK1, PsiLeft, QShift, Op.data(), PsiLeft, mask_row(Op, Row))[Row];
         KMatrixPolyType CTop = inject_right_mask(FMatKTop, PsiLeft, QShift, Op.data(), PsiTri, mask_row(Op, Row))[Row];
         for (auto I = CTop.begin(); I != CTop.end(); ++I)
            C[I->first] += ExpIK * I->second;
         KMatrixPolyType CBot = inject_right_mask(FMatKBot, PsiTri, QShift, Op.data(), PsiLeft, mask_row(Op, Row))[Row];
         for (auto I = CBot.begin(); I != CBot.end(); ++I)
            C[I->first] += std::conj(ExpIK) * I->second;
         KMatrixPolyType CInit = inject_right_mask(FMatK0, PsiTri, QShift, Op.data(), PsiTri, mask_row(Op, Row))[Row];
         for (auto I = CInit.begin(); I != CInit.end(); ++I)
            C[I->first] += I->second;
      }

      // Now do the classification, based on the properties of the diagonal operator.
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;

         FMatK1[Row] = SolveZeroDiagonal(C);
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
         if (Classification.is_unitary() || Classification.is_unclassified())
         {
            // Find the eigenvectors of the transfer matrix: if the operator is
            // proportional to the identity in the scalar sector, then the
            // eigenvectors are those of the usual transfer matrix (TLeft and
            // TRight), otherwise, we must solve for the generalised transfer
            // matrix eigenvalues.
            if (!is_scalar(Diag.Basis1()[0]) || !(Classification.is_complex_identity() || Classification.is_unclassified()))
            {
               if (Verbose > 0)
                  std::cerr << "Solving unitary diagonal component" << std::endl;

               // Find the largest eigenvalue.
               std::complex<double> EValue;
               if (Mode == "initial")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiRight, PsiRight, QShift, ProductMPO(Diag), Tol, Verbose);
               else if (Mode == "top" || Mode == "top-ea")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiLeft, PsiRight, QShift, ProductMPO(Diag), Tol, Verbose);
               else if (Mode == "bottom")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiRight, PsiLeft, QShift, ProductMPO(Diag), Tol, Verbose);
               else if (Mode == "final")
                  std::tie(EValue, TransferEVLeft, TransferEVRight) = get_transfer_eigenpair(PsiLeft, PsiLeft, QShift, ProductMPO(Diag), Tol, Verbose);

               TransferEVRight = delta_shift(TransferEVRight, QShift);

               //EValue = std::conj(EValue); // left eigenvalue, so conjugate (see comment at operator_actions.h)

               if (Verbose > 0)
                  std::cerr << "Eigenvalue of unitary operator is " << EValue << std::endl;

               // If the magnitude of the largest eigenvalue is less than one, we don't need the eigenvectors.
               if (std::abs(norm_frob(Factor) - 1.0) < UnityEpsilon)
               {
                  TransferEVLeft = MatrixOperator();
                  TransferEVRight = MatrixOperator();
               }
               else
                  Factor = EValue;
            }
            else if (Verbose > 0 && Classification.is_identity())
               std::cerr << "Diagonal component is the identity" << std::endl;
            else if (Verbose > 0 && Classification.is_complex_identity())
               std::cerr << "Diagonal component is proportional to the identity" << std::endl;
            else if (Classification.is_unclassified())
            {
               if (Verbose > 0)
                  std::cerr << "Diagonal component is unclassified" << std::endl;
               Factor = 1.0;
            }

            // Multiply by the EA momentum factors for mixed transfer matrices.
            if (Mode == "top" || Mode == "top-ea")
               Factor *= ExpIK;
            else if (Mode == "bottom")
               Factor *= std::conj(ExpIK);

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

            if (Mode == "initial")
               F = DecomposePerpendicularPartsRight(C, Diag, TransferEVRight, TransferEVLeft,
                                                    PsiRight, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "top-ea" && Row == 0)
               // For the EA algorithm, we only need the zero momentum components for the first row.
               F[1.0] = DecomposePerpendicularPartsRight(C[1.0], 1.0, std::conj(ExpIK)*Diag, TransferEVRight, TransferEVLeft,
                                                         PsiLeft, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "top" || Mode == "top-ea")
               F = DecomposePerpendicularPartsRight(C, std::conj(ExpIK)*Diag, TransferEVRight, TransferEVLeft,
                                                    PsiLeft, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "bottom")
               F = DecomposePerpendicularPartsRight(C, ExpIK*Diag, TransferEVRight, TransferEVLeft,
                                                    PsiRight, PsiLeft, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else if (Mode == "final")
               F = DecomposePerpendicularPartsRight(C, Diag, TransferEVRight, TransferEVLeft,
                                                    PsiLeft, PsiLeft, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
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
         FMatK1[Row] = F;
      }
   }
}

void
SolveHamiltonianMPO_EA_Left(StateComponent& E1, std::vector<KMatrixPolyType> const& EMatK0,
                            LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                            LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                            BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                            MatrixOperator const& TRight, std::complex<double> ExpIK,
                            double Tol, int Verbose)
{
   if (E1.is_null())
      E1 = StateComponent(Op.Basis(), PsiRight.Basis1(), PsiLeft.Basis1());
   std::vector<KMatrixPolyType> EMatK1(E1.size());
   for (int i = 0; i < E1.size(); ++i)
   {
      if (!E1[i].is_null())
         EMatK1[i][1.0][0] = E1[i];
   }
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   SolveMPO_EA_Left(EMatK1, EMatK0, std::vector<KMatrixPolyType>(), std::vector<KMatrixPolyType>(),
                    PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK,
                    true, 0, Tol, UnityEpsilon, "top-ea", Verbose);
   for (int i = 0; i < E1.size(); ++i)
   {
      E1[i] = EMatK1[i][1.0].coefficient(0);
   }
   // TODO: Check for diverging components?
}

void
SolveHamiltonianMPO_EA_Right(StateComponent& F1, std::vector<KMatrixPolyType> const& FMatK0,
                             LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                             LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                             BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                             MatrixOperator const& TRight, std::complex<double> ExpIK,
                             double Tol, int Verbose)
{
   if (F1.is_null())
      F1 = StateComponent(Op.Basis(), PsiLeft.Basis2(), PsiRight.Basis2());
   std::vector<KMatrixPolyType> FMatK1(F1.size());
   for (int i = 0; i < F1.size(); ++i)
   {
      if (!F1[i].is_null())
         FMatK1[i][1.0][0] = F1[i];
   }
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   SolveMPO_EA_Right(FMatK1, FMatK0, std::vector<KMatrixPolyType>(), std::vector<KMatrixPolyType>(),
                     PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK,
                     true, 0, Tol, UnityEpsilon, "top-ea", Verbose);
   for (int i = 0; i < F1.size(); ++i)
   {
      F1[i] = FMatK1[i][1.0].coefficient(0);
   }
   // TODO: Check for diverging components?
}

void
SolveHamiltonianMPO_EA_Left(StateComponent& E1, StateComponent const& E0,
                            LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                            LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                            BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                            MatrixOperator const& TRight, std::complex<double> ExpIK,
                            double Tol, int Verbose)
{
   std::vector<KMatrixPolyType> EMatK0(E0.size());
   for (int i = 0; i < E0.size(); ++i)
   {
      if (!E0[i].is_null())
         EMatK0[i][1.0][0] = E0[i];
   }
   SolveHamiltonianMPO_EA_Left(E1, EMatK0, PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK, Tol, Verbose);
}

void
SolveHamiltonianMPO_EA_Right(StateComponent& F1, StateComponent const& F0,
                             LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                             LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                             BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                             MatrixOperator const& TRight, std::complex<double> ExpIK,
                             double Tol, int Verbose)
{
   std::vector<KMatrixPolyType> FMatK0(F0.size());
   for (int i = 0; i < F0.size(); ++i)
   {
      if (!F0[i].is_null())
         FMatK0[i][1.0][0] = F0[i];
   }
   SolveHamiltonianMPO_EA_Right(F1, FMatK0, PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK, Tol, Verbose);
}
