/ -*- C++ -*-
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

#include "triangular_mpo_solver.h"
#include "triangular_mpo_solver_helpers.h"

std::complex<double>
SolveFirstOrderMPO_Left(StateComponent& E, LinearWavefunction const& Psi,
                        QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                        MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (E.is_null())
      E = Initial_E(Op, Psi.Basis1());
   CHECK_EQUAL(E.Basis1(), Psi.Basis1());
   CHECK_EQUAL(E.Basis2(), Psi.Basis1());
   CHECK_EQUAL(E.LocalBasis(), Op.Basis1());
   CHECK_EQUAL(Rho.Basis1(), Psi.Basis1());
   CHECK_EQUAL(Rho.Basis2(), Psi.Basis1());

   MatrixOperator Ident = E[0];

   // The UnityEpsilon is just a paranoid check here, as we don't support
   // eigenvalue 1 on the diagonal
   double UnityEpsilon = 1E-12;

   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Left: fatal: MPO(0,0) must be the identity operator.\n";
      PANIC("Fatal");
   }

   int Dim = Op.Basis1().size();       // dimension of the MPO
   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Left: dimension is " << Dim << std::endl;

   // Column 0 (E[0]) is the Identity, we don't need to solve for it
   for (int Col = 1; Col < Dim-1; ++Col)
   {
      std::vector<std::vector<int> > Mask = mask_column(Op, Col);
      MatrixOperator C = inject_left_mask(E, Psi, Op.data(), Psi, Mask)[Col];
      C = delta_shift(C, QShift);

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << (Col+1) << std::endl;
         E[Col] = C;
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << (Col+1) << std::endl;

         // non-zero diagonal element.  The only case that we support here is
         // an operator with spectral radius strictly < 1
         // if (Classification.is_unitary())
         // {
         //    std::cerr << "SolveSimpleMPO_Left: Unitary operator on the diagonal is not supported!\n";
         //    PANIC("Fatal: unitary")(Col);
         // }

         // Initial guess for linear solver
         if (E[Col].is_null())
            E[Col] = C;

         //    if (Col == 2)
         //       TRACE(C);

         LinearSolve(E[Col], OneMinusTransferLeft_Ortho(Psi, QShift, Diag, Psi, Ident, Rho, false), C, Tol, Verbose);
         //LinearSolve(E[Col], OneMinusTransferLeft(Diag, Psi, QShift), C, Ident, Rho, Tol, Verbose);
      }
   }
   // Final column, must be identity
   int const Col = Dim-1;
   if (!classify(Op(Col,Col), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Left: fatal: MPO(d,d) must be the identity operator.\n";
      PANIC("Fatal");
   }

   std::vector<std::vector<int> > Mask = mask_column(Op, Col);
   MatrixOperator C = inject_left_mask(E, Psi, Op.data(), Psi, Mask)[Col];
   C = delta_shift(C, QShift);

   // The component in the direction of the identity is proportional to the energy
   std::complex<double> Energy = inner_prod(Rho, C);

   // orthogonalize
   C -= Energy * Ident;

   // solve for the first component
   SubProductLeftProject ProdL(Psi, QShift, Rho, Ident);
   if (E[Col].is_null())
      E[Col] = C;
   LinearSolve(E[Col], ProdL, C, Tol, Verbose);

   // Make it Hermitian
   //   E.back() = 0.5 * (E.back() + adjoint(E.back()));

   DEBUG_CHECK(norm_frob(inner_prod(Rho, E.back())) < 1E-12); // make sure that the final operator is orthogonal to the identity

#if !defined(NDEBUG)
   {
      std::vector<KMatrixPolyType> CheckEMat;
      SolveMPO_Left(CheckEMat, Psi, QShift, Op, Ident, Rho, true);
      ComplexPolyType EValues = ExtractOverlap(CheckEMat.back()[1.0], Rho);
      TRACE(EValues);
      MatrixOperator HCheck = CheckEMat.back()[1.0][0];
      TRACE("H matrix elements check")(inner_prod(HCheck - E.back(), Rho));
      TRACE(inner_prod(HCheck, Rho));
   }
#endif

   return Energy;
}

std::complex<double>
SolveFirstOrderMPO_Left(StateComponent& E, InfiniteWavefunctionLeft const& Psi,
                        BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(PsiLinear, Lambda) = get_left_canonical(Psi);
   MatrixOperator Rho = Lambda*Lambda;

   return SolveFirstOrderMPO_Left(E, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}

//
// SolveSimpleMPO_Right
//

std::complex<double>
SolveFirstOrderMPO_Right(StateComponent& F, LinearWavefunction const& Psi,
                         QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                         MatrixOperator const& Rho, double Tol, int Verbose)
{
   if (F.is_null())
      F = Initial_F(Op, Psi.Basis2());
   CHECK_EQUAL(F.Basis1(), Psi.Basis2());
   CHECK_EQUAL(F.Basis2(), Psi.Basis2());
   CHECK_EQUAL(F.LocalBasis(), Op.Basis1());
   CHECK_EQUAL(Rho.Basis1(), Psi.Basis2());
   CHECK_EQUAL(Rho.Basis2(), Psi.Basis2());

   MatrixOperator Ident = F.back();

   // The UnityEpsilon is just a paranoid check here, as we don't support
   // eigenvalue 1 on the diagonal
   double UnityEpsilon = 1E-12;

   if (!classify(Op(0,0), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveFirstOrderMPO_Right: fatal: MPO(0,0) must be the identity operator.\n";
      PANIC("Fatal");
   }

   int Dim = Op.Basis1().size();       // dimension of the MPO
   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Right: dimension is " << Dim << std::endl;

   // Row Dim-1 (F[Dim-1]) is the Identity, we don't need to solve for it
   for (int Row = Dim-2; Row >= 1; --Row)
   {
      std::vector<std::vector<int> > Mask = mask_row(Op, Row);
      MatrixOperator C = inject_right_mask(F, Psi, Op.data(), Psi, Mask)[Row];
      C.delta_shift(adjoint(QShift));

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Row)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;
         F[Row] = C;
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at row " << Row << std::endl;

         // non-zero diagonal element.  The only case that we support here is
         // an operator with spectral radius strictly < 1
         // if (Classification.is_unitary())
         // {
         //    std::cerr << "SolveSimpleMPO_Right: Unitary operator on the diagonal is not supported!\n";
         //    PANIC("Fatal: unitary")(Row);
         // }

         // Initial guess for linear solver
         if (F[Row].is_null())
            F[Row] = C;

         LinearSolve(F[Row], OneMinusTransferRight(Diag, Psi, QShift), C, Tol, Verbose);
      }
   }
   // Final row, must be identity
   int const Row = 0;
   if (!classify(Op(Row,Row), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Right: fatal: MPO(d,d) must be the identity operator.\n";
      PANIC("Fatal");
   }

   std::vector<std::vector<int> > Mask = mask_row(Op, Row);
   MatrixOperator C = inject_right_mask(F, Psi, Op.data(), Psi, Mask)[Row];
   C.delta_shift(adjoint(QShift));

   // The component in the direction of the identity is proportional to the energy
   std::complex<double> Energy = inner_prod(Rho, C);
   // orthogonalize
   C -= Energy * Ident;

   // solve for the first component
   SubProductRightProject ProdR(Psi, QShift, Rho, Ident);
   if (F[Row].is_null())
      F[Row] = C;
   LinearSolve(F[Row], ProdR, C, Tol, Verbose);

   // Make it Hermitian
   //   F.front() = 0.5 * (F.front() + adjoint(F.front()));

   // stability fix
   if (Verbose > 0)
      std::cerr << "Overall constant " << inner_prod(F.front(), F.back()) << '\n';
   F.front() -= inner_prod(F.back(), F.front()) * F.back();

   // remove the spurious constant term from the energy
   if (Verbose > 0)
      std::cerr << "Spurius constant " << inner_prod(F.front(), Rho) << '\n';
   F.front() -= inner_prod(Rho, F.front()) * F.back();

   // Everything here is in the Hermitian representation, so the actual energy is
   // the conjugate
   return std::conj(Energy);
}

std::complex<double>
SolveFirstOrderMPO_Right(StateComponent& F, InfiniteWavefunctionRight const& Psi,
                         BasicTriangularMPO const& Op, double Tol, int Verbose)
{
   LinearWavefunction PsiLinear;
   RealDiagonalOperator Lambda;
   std::tie(Lambda, PsiLinear) = get_right_canonical(Psi);
   MatrixOperator Rho = Lambda*Lambda;
   return SolveSimpleMPO_Right(F, PsiLinear, Psi.qshift(), Op, Rho, Tol, Verbose);
}
