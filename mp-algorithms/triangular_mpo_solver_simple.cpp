// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/triangular_mpo_solver_simple.cpp
//
// Copyright (C) 2009-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

//
// MPO solver for 'simple' triangular MPO's.  These are a bit mis-named,
// they are not so simple, but refers to an MPO that has no momentum dependence,
// and if there are any diagonal components they are proportional to the identity,
// with either prefactor 1.0 or |prefactor| < 1.0

#include "triangular_mpo_solver.h"
#include "triangular_mpo_solver_helpers.h"

void
SolveSimpleMPO_Left(std::vector<MatrixPolyType>& EMat,
                    LinearWavefunction const& Psi, QuantumNumber const& QShift,
                    BasicTriangularMPO const& Op,
                    MatrixOperator const& Identity,
                    MatrixOperator const& Rho, bool NeedFinalMatrix,
                    int Degree, double Tol,
                    double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(Rho.Basis1(), Psi.Basis1());
   CHECK_EQUAL(Rho.Basis2(), Psi.Basis1());
   CHECK_EQUAL(Identity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(Identity.Basis2(), Psi.Basis1());

   //DEBUG_TRACE(Verbose)(Degree)(Tol);

   int Dim = Op.Basis1().size();       // dimension of the MPO
   EMat.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Left: dimension is " << Dim << std::endl;

   // Make sure the (0,0) part is identity
   //DEBUG_TRACE(UnityEpsilon);
   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      // the (0,0) component isn't the identity operator, which is a fatal error.
      // Show some diagnosics and quit.
      std::cerr << "SolveSimpleMPO_Left: fatal: (0,0) component of the MPO must be the identity operator.\n";
      PANIC("Fatal");
   }

   if (EMat[0].empty())
   {
      // Initialize the first E matrix.  These are operators acting in the Basis1()
      EMat[0] = MatrixPolyType(Identity);
   }

   int StartCol = 1;

   // solve recursively column 1 onwards
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving column " << (Col) << " of [0:" << (Dim-1) << "]\n";
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      MatrixPolyType C = inject_left_mask(EMat, Psi, QShift, Op.data(), Psi, mask_column(Op, Col))[Col];

      TRACE(C);

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         //DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << Col << std::endl;
         EMat[Col] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << Col << std::endl;

         ComplexPolyType EParallel;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            if (!Classification.is_identity())
            {
               std::cerr << "SolveSimpleMPO_Left: fatal: unitary (non-identity) component on the diagonal is not supported\n";
               abort();
            }
            HasEigenvalue1 = true;
            //DEBUG_TRACE(UnitMatrixLeft)(UnitMatrixRight);
            if (Verbose > 0)
               std::cerr << "Decomposing parts parallel to the unit matrix\n";
            EParallel = DecomposeParallelParts(C, Identity, Rho, UnityEpsilon, Degree);
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1

         MatrixPolyType E;
         // if we are on the last column and we don't need the matrix elements, then we can
         // skip this operation
         if (Col < Dim-1 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix\n";
            E = DecomposePerpendicularPartsLeft(C, 1.0, Diag, Identity, Rho,
                                                Psi, Psi, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
         {
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last column.\n";
         }

         // Reinsert the components parallel to the unit matrix (if any)
         for (ComplexPolyType::const_iterator J = EParallel.begin(); J != EParallel.end(); ++J)
         {
            // Conj here because this comes from an overlap(x, RightUnitMatrix)
            E[J->first] += std::conj(J->second) * Identity;
         }

         // Finally, set the E matrix element at this column
         //DEBUG_TRACE(E[1.0]);
         EMat[Col] = E;
      }
   }
}

void
SolveSimpleMPO_Right(std::vector<MatrixPolyType>& FMat,
                     LinearWavefunction const& Psi, QuantumNumber const& QShift,
                     BasicTriangularMPO const& Op,
                     MatrixOperator const& Identity,
                     MatrixOperator const& Rho, bool NeedFinalMatrix,
                     int Degree, double Tol,
                     double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(Identity.Basis1(), Psi.Basis2());
   CHECK_EQUAL(Identity.Basis2(), Psi.Basis2());
   CHECK_EQUAL(Rho.Basis1(), Psi.Basis2());
   CHECK_EQUAL(Rho.Basis2(), Psi.Basis2());

   //DEBUG_TRACE(Verbose)(Degree)(Tol);

   int Dim = Op.Basis1().size();       // dimension of the MPO
   FMat.resize(Dim);
   int Row = Dim-1;

   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Right: dimension is " << Dim << std::endl;

   if (!classify(Op(Row,Row), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Right: fatal: MPO(last,last) must be the identity operator.\n";
      PANIC("Fatal");
   }

   // Initialize the first F matrix.  These are operators acting in the Basis1()
   if (FMat[Row].empty())
      FMat[Row] = MatrixPolyType(Identity);

   for (Row = Dim-2; Row >= 0; --Row)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving row " << (Row) << " of [0:" << (Dim-1) << "]\n";
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      MatrixPolyType C = inject_right_mask(FMat, Psi, QShift, Op.data(), Psi, mask_row(Op, Row))[Row];
      //MatrixPolyType C = inject_right(FMat, Psi, QShift, Op.data(), Psi)[Row];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         //DEBUG_TRACE("Zero diagonal element")(Row)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;
         FMat[Row] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at row " << Row << std::endl;

         ComplexPolyType FParallel;

         bool HasEigenvalue1 = false;
         if (Classification.is_unitary())
         {
            if (!Classification.is_identity())
            {
               std::cerr << "SolveSimpleMPO_Right: fatal: unitary (non-identity) component on the diagonal is not supported\n";
               abort();
            }
            HasEigenvalue1 = true;

            if (Verbose > 0)
               std::cerr << "Decomposing parts parallel to the unit matrix\n";
            FParallel = DecomposeParallelParts(C, Identity, Rho, UnityEpsilon, Degree);
         }

         // Now the remaining components, which is anything that is not proportional
         // to an eigenvector of magnitude 1

         MatrixPolyType F;
         // if we are on the last row and we don't need the matrix elements, then we can
         // skip this operation
         if (Row > 0 || NeedFinalMatrix)
         {
            if (Verbose > 0)
               std::cerr << "Decomposing parts perpendicular to the unit matrix\n";
            F = DecomposePerpendicularPartsRight(C, 1.0, Diag, Identity, Rho,
                                                 Psi, Psi, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
         }
         else if (Verbose > 0)
         {
            std::cerr << "Skipping parts perpendicular to the unit matrix for the last row.\n";
         }

         // Reinsert the components parallel to the unit matrix (if any)
         for (ComplexPolyType::const_iterator J = FParallel.begin(); J != FParallel.end(); ++J)
         {
            // Conj here because this comes from an overlap(x, LeftUnitMatrix)
            F[J->first] += std::conj(J->second) * Identity;
         }

         // Finally, set the E matrix element at this column
         //DEBUG_TRACE(E[1.0]);
         FMat[Row] = F;
      }
   }
}

std::complex<double>
SolveHamiltonianMPO_Left(StateComponent& E, LinearWavefunction const& Psi,
   QuantumNumber const& QShift, BasicTriangularMPO const& Op,
   MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (E.is_null())
      E = Initial_E(Op, Psi.Basis1());
   std::vector<MatrixPolyType> EMat(E.size());
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
      //DEBUG_TRACE(Remainder);
      //DEBUG_TRACE(norm_frob(Remainder));
      //DEBUG_TRACE(inner_prod(Rho, Remainder));
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
SolveHamiltonianMPO_Left(StateComponent& E, InfiniteWavefunctionLeft const& Psi,
   BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(PsiLinear, Lambda) = get_left_canonical(Psi);
   MatrixOperator Rho = delta_shift(Lambda*Lambda, Psi.qshift());
   return SolveHamiltonianMPO_Left(E, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}

std::complex<double>
SolveHamiltonianMPO_Right(StateComponent& F, LinearWavefunction const& Psi,
   QuantumNumber const& QShift, BasicTriangularMPO const& Op,
   MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (F.is_null())
      F = Initial_F(Op, Psi.Basis1());
   std::vector<MatrixPolyType> FMat(F.size());
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
      //DEBUG_TRACE(Remainder);
      //DEBUG_TRACE(norm_frob(Remainder));
      //DEBUG_TRACE(inner_prod(Rho, Remainder));
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
SolveHamiltonianMPO_Right(StateComponent& F, InfiniteWavefunctionRight const& Psi,
   BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(Lambda, PsiLinear) = get_right_canonical(Psi);
   MatrixOperator Rho = delta_shift(Lambda*Lambda, adjoint(Psi.qshift()));
   return SolveHamiltonianMPO_Right(F, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}
