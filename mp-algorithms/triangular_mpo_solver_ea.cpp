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

#include "triangular_mpo_solver.h"
#include "triangular_mpo_solver_helpers.h"

void
SolveSimpleMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK1, std::vector<MatrixPolyType> const& EMat0,
                       LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                       LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                       BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                       MatrixOperator const& TRight, std::complex<double> ExpIK,
                       bool NeedFinalMatrix, int Degree, double Tol,
                       double UnityEpsilon, int Verbose, bool Bottom)
{
   DEBUG_TRACE(Verbose)(Degree)(Tol);

   if (Bottom)
      ExpIK = std::conj(ExpIK);

   int Dim = Op.Basis1().size();       // dimension of the MPO
   EMatK1.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_EA_Left: dimension is " << Dim << std::endl;

   // Make sure the (0,0) part is identity
   DEBUG_TRACE(UnityEpsilon);
   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      // the (0,0) component isn't the identity operator, which is a fatal error.
      // Show some diagnosics and quit.
      std::cerr << "SolveSimpleMPO_EA_Left: fatal: (0,0) component of the MPO must be the identity operator.\n";
      PANIC("Fatal");
   }

   int StartCol = 0;

   // solve recursively column 0 onwards
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving column " << (Col) << " of [0:" << (Dim-1) << "]\n";
      }

      // Generate the next C matrices
      KMatrixPolyType C;
      if (!Bottom)
         C = inject_left_mask(EMatK1, PsiRight, QShift, Op.data(), PsiLeft, mask_column(Op, Col))[Col];
      else
         C = inject_left_mask(EMatK1, PsiLeft, QShift, Op.data(), PsiRight, mask_column(Op, Col))[Col];
      for (auto I = C.begin(); I != C.end(); ++I)
         I->second *= ExpIK;
      if (!Bottom)
         C[1.0] += inject_left_mask(EMat0, PsiTri, QShift, Op.data(), PsiLeft, mask_column(Op, Col))[Col];
      else
         C[1.0] += inject_left_mask(EMat0, PsiLeft, QShift, Op.data(), PsiTri, mask_column(Op, Col))[Col];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         EMatK1[Col] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

         KComplexPolyType EParallel;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            if (!Classification.is_identity())
            {
               std::cerr << "SolveSimpleMPO_EA_Left: fatal: unitary (non-identity) component on the diagonal is not supported\n";
               abort();
            }
            // We only need to solve for parallel component if PsiLeft and
            // PsiRight are the same state.
            if (!TLeft.is_null())
            {
               HasEigenvalue1 = true;
               //DEBUG_TRACE(UnitMatrixLeft)(UnitMatrixRight);
               if (Verbose > 0)
                  std::cerr << "Decomposing parts parallel to the unit matrix\n";
               EParallel = DecomposeParallelPartsWithMomentum(C, ExpIK, TLeft, TRight, UnityEpsilon, Degree);
            }
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1

         KMatrixPolyType E;
         // if we are on the last column and we don't need the matrix elements, then we can
         // skip this operation
         if (Col < Dim-1 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix\n";
            if (!Bottom)
               E = DecomposePerpendicularPartsLeft(C, ExpIK*Diag, TLeft, TRight, PsiRight, PsiLeft, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
            else
               E = DecomposePerpendicularPartsLeft(C, ExpIK*Diag, TLeft, TRight, PsiLeft, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
         {
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last column.\n";
         }

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               // Conj here because this comes from an overlap(x, RightUnitMatrix)
               E[I->first][J->first] += std::conj(J->second) * TLeft;
            }
         }

         // Finally, set the E matrix element at this column
         EMatK1[Col] = E;
      }
   }
}

void
SolveSimpleMPO_EA_Left_Final(std::vector<KMatrixPolyType>& EMatK1, std::vector<MatrixPolyType> const& EMat0,
                       std::vector<KMatrixPolyType> const& EMatKTop, std::vector<KMatrixPolyType> const& EMatKBot,
                       LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                       LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                       BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                       MatrixOperator const& TRight, std::complex<double> ExpIK,
                       bool NeedFinalMatrix, int Degree, double Tol,
                       double UnityEpsilon, int Verbose)
{
   DEBUG_TRACE(Verbose)(Degree)(Tol);

   int Dim = Op.Basis1().size();       // dimension of the MPO
   EMatK1.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_EA_Left: dimension is " << Dim << std::endl;

   // Make sure the (0,0) part is identity
   DEBUG_TRACE(UnityEpsilon);
   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      // the (0,0) component isn't the identity operator, which is a fatal error.
      // Show some diagnosics and quit.
      std::cerr << "SolveSimpleMPO_EA_Left: fatal: (0,0) component of the MPO must be the identity operator.\n";
      PANIC("Fatal");
   }

   int StartCol = 0;

   // solve recursively column 0 onwards
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving column " << (Col) << " of [0:" << (Dim-1) << "]\n";
      }

      // Generate the next C matrices
      KMatrixPolyType C = inject_left_mask(EMatK1, PsiRight, QShift, Op.data(), PsiRight, mask_column(Op, Col))[Col];
      KMatrixPolyType CTop = inject_left_mask(EMatKTop, PsiRight, QShift, Op.data(), PsiTri, mask_column(Op, Col))[Col];
      for (auto I = CTop.begin(); I != CTop.end(); ++I)
         C[I->first] += ExpIK * I->second;
      KMatrixPolyType CBot = inject_left_mask(EMatKBot, PsiTri, QShift, Op.data(), PsiRight, mask_column(Op, Col))[Col];
      for (auto I = CBot.begin(); I != CBot.end(); ++I)
         C[I->first] += std::conj(ExpIK) * I->second;
      C[1.0] += inject_left_mask(EMat0, PsiTri, QShift, Op.data(), PsiTri, mask_column(Op, Col))[Col];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         EMatK1[Col] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

         KComplexPolyType EParallel;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            if (!Classification.is_identity())
            {
               std::cerr << "SolveSimpleMPO_EA_Left: fatal: unitary (non-identity) component on the diagonal is not supported\n";
               abort();
            }
            HasEigenvalue1 = true;
            //DEBUG_TRACE(UnitMatrixLeft)(UnitMatrixRight);
            if (Verbose > 0)
               std::cerr << "Decomposing parts parallel to the unit matrix\n";
            EParallel = DecomposeParallelPartsWithMomentum(C, 1.0, TLeft, TRight, UnityEpsilon, Degree);
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1

         KMatrixPolyType E;
         // if we are on the last column and we don't need the matrix elements, then we can
         // skip this operation
         if (Col < Dim-1 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix\n";
            E = DecomposePerpendicularPartsLeft(C, Diag, TLeft, TRight, PsiRight, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
         {
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last column.\n";
         }

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               // Conj here because this comes from an overlap(x, RightUnitMatrix)
               E[I->first][J->first] += std::conj(J->second) * TLeft;
            }
         }

         // Finally, set the E matrix element at this column
         EMatK1[Col] = E;
      }
   }
}

void
SolveSimpleMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK1, std::vector<MatrixPolyType> const& FMat0,
                        LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                        LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                        BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                        MatrixOperator const& TRight, std::complex<double> ExpIK,
                        bool NeedFinalMatrix, int Degree, double Tol,
                        double UnityEpsilon, int Verbose)
{
   DEBUG_TRACE(Verbose)(Degree)(Tol);

   int Dim = Op.Basis1().size();       // dimension of the MPO
   FMatK1.resize(Dim);
   int Row = Dim-1;

   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_EA_Right: dimension is " << Dim << std::endl;

   if (!classify(Op(Row,Row), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_EA_Right: fatal: MPO(last,last) must be the identity operator.\n";
      PANIC("Fatal");
   }

   for (Row = Dim-1; Row >= 0; --Row)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving row " << (Row) << " of [0:" << (Dim-1) << "]\n";
      }

      // Generate the next C matrices
      KMatrixPolyType C = inject_right_mask(FMatK1, PsiLeft, QShift, Op.data(), PsiRight, mask_row(Op, Row))[Row];
      for (auto I = C.begin(); I != C.end(); ++I)
         I->second *= ExpIK;
      C[1.0] += inject_right_mask(FMat0, PsiTri, QShift, Op.data(), PsiRight, mask_row(Op, Row))[Row];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Row)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;
         FMatK1[Row] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at row " << Row << std::endl;

         KComplexPolyType FParallel;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            if (!Classification.is_identity())
            {
               std::cerr << "SolveSimpleMPO_EA_Right: fatal: unitary (non-identity) component on the diagonal is not supported\n";
               abort();
            }
            // We only need to solve for parallel component if PsiLeft and
            // PsiRight are the same state.
            if (!TLeft.is_null())
            {
               HasEigenvalue1 = true;

               if (Verbose > 0)
                  std::cerr << "Decomposing parts parallel to the unit matrix\n";
               FParallel = DecomposeParallelPartsWithMomentum(C, ExpIK, TRight, TLeft, UnityEpsilon, Degree);
            }
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1

         KMatrixPolyType F;
         // if we are on the last row and we don't need the matrix elements, then we can
         // skip this operation
         if (Row > 0 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix\n";
            F = DecomposePerpendicularPartsRight(C, std::conj(ExpIK)*Diag, TRight, TLeft, PsiLeft, PsiRight, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
         {
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last row.\n";
         }

         // Reinsert the components parallel to the unit matrix (if any)
         for (KComplexPolyType::const_iterator I = FParallel.begin(); I != FParallel.end(); ++I)
         {
            for (ComplexPolyType::const_iterator J = I->second.begin(); J != I->second.end(); ++J)
            {
               // Conj here because this comes from an overlap(x, LeftUnitMatrix)
               F[I->first][J->first] += std::conj(J->second) * TRight;
            }
         }

         // Finally, set the F matrix element at this row
         FMatK1[Row] = F;
      }
   }
}

std::complex<double>
SolveHamiltonianMPO_Left(std::vector<MatrixPolyType>& EMat, StateComponent& E, LinearWavefunction const& Psi,
                         QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                         MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (E.is_null())
      E = Initial_E(Op, Psi.Basis1());
   EMat = std::vector<MatrixPolyType>(E.size());
   for (int i = 0; i < E.size(); ++i)
   {
      if (!E[i].is_null())
         EMat[i][0] = E[i];
   }
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   SolveSimpleMPO_Left(EMat, Psi, QShift, Op, E.front(), Rho, true, 0, Tol, UnityEpsilon, Verbose);
   for (int i = 0; i < E.size(); ++i)
   {
      E[i] = EMat[i].coefficient(0);
   }
   std::complex<double> Energy = inner_prod(Rho, EMat.back()[1]);
   // Check that the linear part of the Hamiltonian is a constant
   MatrixOperator Remainder = EMat.back()[1] - Energy*E.front();
   if (norm_frob(Remainder) > Tol * norm_frob(Energy))
   {
      std::cerr << "SolveHamiltonianMPO_Left: warning: Hamiltonian has diverging matrix elements.\n";
      std::cerr << "Norm of remainder = " << norm_frob(Remainder) << '\n';
      DEBUG_TRACE(Remainder);
      DEBUG_TRACE(norm_frob(Remainder));
      DEBUG_TRACE(inner_prod(Rho, Remainder));
      //std::abort();
   }
   if (EMat.back().degree() > 1)
   {
      for (int d = 2; d <= EMat.back().degree(); ++d)
      {
         if (norm_frob(EMat.back().coefficient(d)) > Tol * norm_frob(Energy))
         {
            std::cerr << "SolveHamiltonianMPO_Left: error: energy per site diverges at order " << d << " with component magnitude " << norm_frob(EMat.back().coefficient(d)) << '\n';
            std::abort();
         }
      }
   }
   return Energy;
}

std::complex<double>
SolveHamiltonianMPO_Left(std::vector<MatrixPolyType>& EMat, StateComponent& E, InfiniteWavefunctionLeft const& Psi,
                         BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(PsiLinear, Lambda) = get_left_canonical(Psi);
   MatrixOperator Rho = delta_shift(Lambda*Lambda, Psi.qshift());
   return SolveHamiltonianMPO_Left(EMat, E, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}

std::complex<double>
SolveHamiltonianMPO_Right(std::vector<MatrixPolyType>& FMat, StateComponent& F, LinearWavefunction const& Psi,
                          QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                          MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (F.is_null())
      F = Initial_F(Op, Psi.Basis1());
   FMat = std::vector<MatrixPolyType>(F.size());
   for (int i = 0; i < F.size(); ++i)
   {
      if (!F[i].is_null())
         FMat[i][0] = F[i];
   }
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   SolveSimpleMPO_Right(FMat, Psi, QShift, Op, F.back(), Rho, true, 0, Tol, UnityEpsilon, Verbose);
   for (int i = 0; i < F.size(); ++i)
   {
      F[i] = FMat[i].coefficient(0);
   }
   std::complex<double> Energy = inner_prod(Rho, FMat.front()[1]);
   // Check that the linear part of the Hamiltonian is a constant
   MatrixOperator Remainder = FMat.front()[1] - Energy*F.back();
   if (norm_frob(Remainder) > Tol * norm_frob(Energy))
   {
      std::cerr << "SolveHamiltonianMPO_Right: warning: Hamiltonian has diverging matrix elements.\n";
      std::cerr << "Norm of remainder = " << norm_frob(Remainder) << '\n';
      DEBUG_TRACE(Remainder);
      DEBUG_TRACE(norm_frob(Remainder));
      DEBUG_TRACE(inner_prod(Rho, Remainder));
      //std::abort();
   }
   if (FMat.front().degree() > 1)
   {
      for (int d = 2; d <= FMat.front().degree(); ++d)
      {
         if (norm_frob(FMat.front().coefficient(d)) > Tol * norm_frob(Energy))
         {
            std::cerr << "SolveHamiltonianMPO_Right: error: energy per site diverges at order " << d << " with component magnitude " << norm_frob(FMat.front().coefficient(d)) << '\n';
            TRACE(inner_prod(Rho,FMat.front().coefficient(d)));
            std::abort();
         }
      }
   }
   // Everything here is in the Hermitian representation, so the actual energy is the conjugate
   return std::conj(Energy);
}

std::complex<double>
SolveHamiltonianMPO_Right(std::vector<MatrixPolyType>& FMat, StateComponent& F, InfiniteWavefunctionRight const& Psi,
                          BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(Lambda, PsiLinear) = get_right_canonical(Psi);
   MatrixOperator Rho = delta_shift(Lambda*Lambda, adjoint(Psi.qshift()));
   return SolveHamiltonianMPO_Right(FMat, F, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}

void
SolveHamiltonianMPO_EA_Left(StateComponent& E1, std::vector<MatrixPolyType> const& EMat0,
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
   SolveSimpleMPO_EA_Left(EMatK1, EMat0, PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK,
                          true, 0, Tol, UnityEpsilon, Verbose);
   for (int i = 0; i < E1.size(); ++i)
   {
      E1[i] = EMatK1[i][1.0].coefficient(0);
   }
#if 0  // This will not work if TLeft and TRight are not defined.
   std::complex<double> Energy = inner_prod(TRight, EMat1.back()[1]);
   // Check that the linear part of the Hamiltonian is a constant
   MatrixOperator Remainder = EMat1.back()[1] - Energy*TLeft;
   if (norm_frob(Remainder) > Tol * norm_frob(Energy))
   {
      std::cerr << "SolveHamiltonianMPO_Left: warning: Hamiltonian has diverging matrix elements.\n";
      std::cerr << "Norm of remainder = " << norm_frob(Remainder) << '\n';
      DEBUG_TRACE(Remainder);
      DEBUG_TRACE(norm_frob(Remainder));
      DEBUG_TRACE(inner_prod(TRight, Remainder));
      //std::abort();
   }
   if (EMat1.back().degree() > 1)
   {
      for (int d = 2; d <= EMat1.back().degree(); ++d)
      {
         if (norm_frob(EMat1.back().coefficient(d)) > Tol * norm_frob(Energy))
         {
            std::cerr << "SolveHamiltonianMPO_Left: error: energy per site diverges at order " << d << " with component magnitude " << norm_frob(EMat1.back().coefficient(d)) << '\n';
            std::abort();
         }
      }
   }
#endif
}

void
SolveHamiltonianMPO_EA_Right(StateComponent& F1, std::vector<MatrixPolyType> const& FMat0,
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
   SolveSimpleMPO_EA_Right(FMatK1, FMat0, PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK,
                           true, 0, Tol, UnityEpsilon, Verbose);
   for (int i = 0; i < F1.size(); ++i)
   {
      F1[i] = FMatK1[i][1.0].coefficient(0);
   }
#if 0
   std::complex<double> Energy = inner_prod(TLeft, FMat1.front()[1]);
   // Check that the linear part of the Hamiltonian is a constant
   MatrixOperator Remainder = FMat1.front()[1] - Energy*TRight
   if (norm_frob(Remainder) > Tol * norm_frob(Energy))
   {
      std::cerr << "SolveHamiltonianMPO_Right: warning: Hamiltonian has diverging matrix elements.\n";
      std::cerr << "Norm of remainder = " << norm_frob(Remainder) << '\n';
      DEBUG_TRACE(Remainder);
      DEBUG_TRACE(norm_frob(Remainder));
      DEBUG_TRACE(inner_prod(TLeft, Remainder));
      //std::abort();
   }
   if (FMat1.front().degree() > 1)
   {
      for (int d = 2; d <= FMat1.front().degree(); ++d)
      {
         if (norm_frob(FMat1.front().coefficient(d)) > Tol * norm_frob(Energy))
         {
            std::cerr << "SolveHamiltonianMPO_Right: error: energy per site diverges at order " << d << " with component magnitude " << norm_frob(FMat1.front().coefficient(d)) << '\n';
            TRACE(inner_prod(Rho,FMat1.front().coefficient(d)));
            std::abort();
         }
      }
   }
#endif
}

void
SolveHamiltonianMPO_EA_Left(StateComponent& E1, StateComponent const& E0,
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
   std::vector<MatrixPolyType> EMat0(E0.size());
   for (int i = 0; i < E0.size(); ++i)
   {
      if (!E0[i].is_null())
         EMat0[i][0] = E0[i];
   }
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   SolveSimpleMPO_EA_Left(EMatK1, EMat0, PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK,
                          true, 0, Tol, UnityEpsilon, Verbose);
   for (int i = 0; i < E1.size(); ++i)
   {
      //TRACE(i)(norm_frob(EMatK1[i][1.0].coefficient(0)))(norm_frob(EMatK1[i][ExpIK].coefficient(0)));
      E1[i] = EMatK1[i][1.0].coefficient(0);
   }
}

void
SolveHamiltonianMPO_EA_Right(StateComponent& F1, StateComponent const& F0,
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
   std::vector<MatrixPolyType> FMat0(F0.size());
   for (int i = 0; i < F0.size(); ++i)
   {
      if (!F0[i].is_null())
         FMat0[i][0] = F0[i];
   }
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   SolveSimpleMPO_EA_Right(FMatK1, FMat0, PsiLeft, PsiRight, PsiTri, QShift, Op, TLeft, TRight, ExpIK,
                           true, 0, Tol, UnityEpsilon, Verbose);
   for (int i = 0; i < F1.size(); ++i)
   {
      //TRACE(i)(norm_frob(FMatK1[i][1.0].coefficient(0)))(norm_frob(FMatK1[i][ExpIK].coefficient(0)));
      F1[i] = FMatK1[i][1.0].coefficient(0);
   }
}
