// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/triangular_mpo_solver.cpp
//
// Copyright (C) 2009-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
                    MatrixOperator const& LeftIdentity,
                    MatrixOperator const& RightIdentity, bool NeedFinalMatrix,
                    int Degree, double Tol,
                    double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(RightIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(RightIdentity.Basis2(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis2(), Psi.Basis1());

   DEBUG_TRACE(Verbose)(Degree)(Tol);

   int Dim = Op.Basis1().size();       // dimension of the MPO
   EMat.resize(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_Left: dimension is " << Dim << std::endl;

   // Make sure the (0,0) part is identity
   DEBUG_TRACE(UnityEpsilon);
   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      // the (0,0) component isn't the identity operator, which is a fatal error.
      // Show some diagnosics and quit.
      std::cerr << "SolveSimpleMPO_Left: fatal: (0,0) component of the MPO must be the identity operator.\n";
   }

   if (EMat[0].empty())
   {
      // Initialize the first E matrix.  These are operators acting in the Basis1()
      EMat[0] = MatrixPolyType(LeftIdentity);
   }

   int StartCol = 1;

   // solve recursively column 1 onwards
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving column " << (Col) << " of [0:" << Dim << "]\n";
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      MatrixPolyType C = inject_left_mask(EMat, Psi, QShift, Op.data(), Psi, mask_column(Op, Col))[Col];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
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
            EParallel = DecomposeParallelParts(C, LeftIdentity, RightIdentity, UnityEpsilon, Degree);
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
            E = DecomposePerpendicularPartsLeft(C, 1.0, Diag, LeftIdentity, RightIdentity,
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
            E[J->first] += std::conj(J->second) * LeftIdentity;
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
                     MatrixOperator const& LeftIdentity,
                     MatrixOperator const& RightIdentity, bool NeedFinalMatrix,
                     int Degree, double Tol,
                     double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(RightIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(RightIdentity.Basis2(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis1(), Psi.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis2(), Psi.Basis1());

   DEBUG_TRACE(Verbose)(Degree)(Tol);

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
      FMat[Row] = MatrixPolyType(RightIdentity);

   for (Row = Dim-2; Row >= 0; --Row)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving row " << (Row) << " of [0:" << (Dim-1) << "]\n";
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      MatrixPolyType C = inject_right_mask(FMat, Psi, QShift, Op.data(), Psi, mask_row(Op, Row))[Row];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Row)(Diag);
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
            FParallel = DecomposeParallelParts(C, RightIdentity, LeftIdentity, UnityEpsilon, Degree);
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
            F = DecomposePerpendicularPartsRight(C, 1.0, Diag, LeftIdentity, RightIdentity,
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
            F[J->first] += std::conj(J->second) * RightIdentity;
         }

         // Finally, set the E matrix element at this column
         //DEBUG_TRACE(E[1.0]);
         FMat[Row] = F;
      }
   }
}
