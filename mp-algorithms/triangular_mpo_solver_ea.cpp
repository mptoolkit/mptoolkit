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

// The tolerance for determining whether ExpIK == 1.0.
double const ExpIKTol = 1e-14;

struct OneMinusTransferLeftEA
{
   OneMinusTransferLeftEA(GenericMPO const& Op_, LinearWavefunction const& Psi1_,
                          LinearWavefunction const& Psi2_, QuantumNumber const& QShift_,
                          std::complex<double> ExpIK_)
      : Op(Op_), Psi1(Psi1_), Psi2(Psi2_), QShift(QShift_), ExpIK(ExpIK_)
   {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      return x - ExpIK * delta_shift(inject_left(x, Psi1, Op, Psi2), QShift);
   }

   GenericMPO const& Op;
   LinearWavefunction const& Psi1;
   LinearWavefunction const& Psi2;
   QuantumNumber const& QShift;
   std::complex<double> ExpIK;
};

struct SubProductLeftProjectEA
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeftProjectEA(LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
                           QuantumNumber const& QShift_, MatrixOperator const& Proj_,
                           MatrixOperator const& Ident_, std::complex<double> ExpIK_)
      : Psi1(Psi1_), Psi2(Psi2_), QShift(QShift_), Proj(Proj_), Ident(Ident_), ExpIK(ExpIK_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In; //delta_shift(In, QShift);
      auto I1 = Psi1.begin();
      auto I2 = Psi2.begin();
      while (I1 != Psi1.end())
      {
         Result = operator_prod(herm(*I2), Result, *I1);
         ++I1, ++I2;
      }
      Result *= ExpIK;
      Result = In - delta_shift(Result, QShift);
      //Result = 0.5 * (Result + adjoint(Result));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi1;
   LinearWavefunction const& Psi2;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
   std::complex<double> ExpIK;
};

void
SolveFirstOrderMPO_EA_Left(StateComponent& E1, StateComponent const& E0,
                           LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                           LinearWavefunction const& PsiTri,
                           QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                           MatrixOperator const& Ident, MatrixOperator const& Rho,
                           std::complex<double> ExpIK, double Tol, int Verbose)
{
   if (E1.is_null())
      E1 = StateComponent(Op.Basis(), PsiRight.Basis1(), PsiLeft.Basis1());

   int Dim = Op.Basis1().size();       // dimension of the MPO
   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Left2: dimension is " << Dim << std::endl;

   // First column
   int Col = 0;

   // The UnityEpsilon is just a paranoid check here, as we don't support
   // eigenvalue 1 on the diagonal
   double UnityEpsilon = 1E-12;

   if (!classify(Op(Col,Col), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Left2: fatal: MPO(0,0) must be the identity operator.\n";
      PANIC("Fatal");
   }

   std::vector<std::vector<int> > Mask = mask_column(Op, Col);
   MatrixOperator C = inject_left_mask(E0, PsiTri, Op.data(), PsiLeft, Mask)[Col];
   C = delta_shift(C, QShift);

   // solve for the first component
   // If the spectral radius of the transfer matrix is < 1 or if k != 0, then
   // we do not need to orthogonalize against the leading eigenvector.
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > ExpIKTol)
   {
      E1[Col] = C;
      LinearSolve(E1[Col], OneMinusTransferLeftEA(Op(Col, Col), PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      E1[Col] = C;
      LinearSolve(E1[Col], SubProductLeftProjectEA(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }

   for (Col = 1; Col < Dim-1; ++Col)
   {
      Mask = mask_column(Op, Col);
      C = inject_left_mask(E0, PsiTri, Op.data(), PsiLeft, Mask)[Col]
          + ExpIK * inject_left_mask(E1, PsiRight, Op.data(), PsiLeft, Mask)[Col];
      C = delta_shift(C, QShift);

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << (Col+1) << std::endl;
         E1[Col] = C;
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << (Col+1) << std::endl;

         // non-zero diagonal element.  The only case that we support here is
         // an operator with spectral radius strictly < 1
         if (Classification.is_unitary())
         {
            std::cerr << "SolveSimpleMPO_Left2: Unitary operator on the diagonal is not supported!\n";
            PANIC("Fatal: unitary")(Col);
         }

         // Initial guess for linear solver
         E1[Col] = C;

	 LinearSolve(E1[Col], OneMinusTransferLeftEA(Diag, PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
      }
   }
   // Final column, must be identity
   Col = Dim-1;
   if (!classify(Op(Col,Col), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Left2: fatal: MPO(d,d) must be the identity operator.\n";
      PANIC("Fatal");
   }

   Mask = mask_column(Op, Col);
   C = inject_left_mask(E0, PsiTri, Op.data(), PsiLeft, Mask)[Col]
       + ExpIK * inject_left_mask(E1, PsiRight, Op.data(), PsiLeft, Mask)[Col];
   C = delta_shift(C, QShift);

   // solve for the final component
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > ExpIKTol)
   {
      E1[Col] = C;
      LinearSolve(E1[Col], OneMinusTransferLeftEA(Op(Col, Col), PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      E1[Col] = C;
      LinearSolve(E1[Col], SubProductLeftProjectEA(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}

struct SubProductLeftProjectEAOp
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeftProjectEAOp(GenericMPO const& Op_, LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
                             QuantumNumber const& QShift_, MatrixOperator const& Proj_,
                             MatrixOperator const& Ident_, std::complex<double> ExpIK_)
      : Op(Op_), Psi1(Psi1_), Psi2(Psi2_), QShift(QShift_), Proj(Proj_), Ident(Ident_), ExpIK(ExpIK_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In - ExpIK * delta_shift(inject_left(In, Psi1, Op, Psi2), QShift);
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi1;
   LinearWavefunction const& Psi2;
   GenericMPO const& Op;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
   std::complex<double> ExpIK;
};

void
SolveStringMPO_EA_Left(MatrixOperator& E1, MatrixOperator const& E0,
                       LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                       LinearWavefunction const& PsiTri,
                       QuantumNumber const& QShift, ProductMPO const& Op,
                       MatrixOperator const& Ident, MatrixOperator const& Rho,
                       std::complex<double> ExpIK, double Tol, int Verbose)
{
   CHECK(Op.is_string());

   if (E1.is_null())
      E1 = MatrixOperator(PsiRight.Basis1(), PsiLeft.Basis1());

   MatrixOperator C = inject_left(E0, PsiTri, Op.data(), PsiLeft);
   C = delta_shift(C, QShift);

   if (Rho.is_null() || std::abs(ExpIK - 1.0) > ExpIKTol)
   {
      E1 = C;
      LinearSolve(E1, OneMinusTransferLeftEA(Op, PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      E1 = C;
      LinearSolve(E1, SubProductLeftProjectEAOp(Op, PsiRight, PsiLeft, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}

struct OneMinusTransferRightEA
{
   OneMinusTransferRightEA(GenericMPO const& Op_, LinearWavefunction const& Psi1_,
                           LinearWavefunction const& Psi2_, QuantumNumber const& QShift_,
                           std::complex<double> ExpIK_)
      : Op(Op_), Psi1(Psi1_), Psi2(Psi2_), QShift(QShift_), ExpIK(ExpIK_)
   {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      return x - ExpIK * delta_shift(inject_right(x, Psi1, Op, Psi2), adjoint(QShift));
   }

   GenericMPO const& Op;
   LinearWavefunction const& Psi1;
   LinearWavefunction const& Psi2;
   QuantumNumber const& QShift;
   std::complex<double> ExpIK;
};

struct SubProductRightProjectEA
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRightProjectEA(LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
                            QuantumNumber const& QShift_, MatrixOperator const& Proj_,
                            MatrixOperator const& Ident_, std::complex<double> ExpIK_)
      : Psi1(Psi1_), Psi2(Psi2_), QShift(QShift_), Proj(Proj_), Ident(Ident_), ExpIK(ExpIK_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In;
      auto I1 = Psi1.end();
      auto I2 = Psi2.end();
      while (I1 != Psi1.begin())
      {
         --I1, --I2;
         Result = operator_prod(*I1, Result, herm(*I2));
      }
      Result *= ExpIK;
      Result = delta_shift(Result, adjoint(QShift));
      Result = In - Result;
      //Result = 0.5 * (Result + adjoint(Result));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi1;
   LinearWavefunction const& Psi2;
   QuantumNumber QShift;
   MatrixOperator const& Proj;
   MatrixOperator const& Ident;
   std::complex<double> ExpIK;
};

void
SolveFirstOrderMPO_EA_Right(StateComponent& F1, StateComponent const& F0,
                            LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                            LinearWavefunction const& PsiTri,
                            QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                            MatrixOperator const& Rho, MatrixOperator const& Ident,
                            std::complex<double> ExpIK, double Tol, int Verbose)
{
   if (F1.is_null())
      F1 = StateComponent(Op.Basis(), PsiLeft.Basis2(), PsiRight.Basis2());

   int Dim = Op.Basis1().size();       // dimension of the MPO
   if (Verbose > 0)
      std::cerr << "SolveSimpleMPO_Right2: dimension is " << Dim << std::endl;

   // Final row
   int Row = Dim-1;

   // The UnityEpsilon is just a paranoid check here, as we don't support
   // eigenvalue 1 on the diagonal
   double UnityEpsilon = 1E-12;

   if (!classify(Op(Row,Row), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Right2: fatal: MPO(d,d) must be the identity operator.\n";
      PANIC("Fatal");
   }

   std::vector<std::vector<int> > Mask = mask_row(Op, Row);
   MatrixOperator C = inject_right_mask(F0, PsiTri, Op.data(), PsiRight, Mask)[Row];
   C.delta_shift(adjoint(QShift));

   // solve for the final component
   // If the spectral radius of the transfer matrix is < 1 or if k != 0, then
   // we do not need to orthogonalize against the leading eigenvector.
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > ExpIKTol)
   {
      F1[Row] = C;
      LinearSolve(F1[Row], OneMinusTransferRightEA(Op(Row, Row), PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      F1[Row] = C;
      LinearSolve(F1[Row], SubProductRightProjectEA(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
   
   for (Row = Dim-2; Row >= 1; --Row)
   {
      Mask = mask_row(Op, Row);
      C = inject_right_mask(F0, PsiTri, Op.data(), PsiRight, Mask)[Row]
          + ExpIK * inject_right_mask(F1, PsiLeft, Op.data(), PsiRight, Mask)[Row];
      C.delta_shift(adjoint(QShift));

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Row, Row);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Row)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at row " << Row << std::endl;
         F1[Row] = C;
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at row " << Row << std::endl;

         // non-zero diagonal element.  The only case that we support here is
         // an operator with spectral radius strictly < 1
         if (Classification.is_unitary())
         {
            std::cerr << "SolveSimpleMPO_Right2: Unitary operator on the diagonal is not supported!\n";
            PANIC("Fatal: unitary")(Row);
         }

         // Initial guess for linear solver
         F1[Row] = C;

	 LinearSolve(F1[Row], OneMinusTransferRightEA(Diag, PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
      }
   }
   // First row, must be identity
   Row = 0;
   if (!classify(Op(Row,Row), UnityEpsilon).is_identity())
   {
      std::cerr << "SolveSimpleMPO_Right2: fatal: MPO(0,0) must be the identity operator.\n";
      PANIC("Fatal");
   }

   Mask = mask_row(Op, Row);
   C = inject_right_mask(F0, PsiTri, Op.data(), PsiRight, Mask)[Row]
       + ExpIK * inject_right_mask(F1, PsiLeft, Op.data(), PsiRight, Mask)[Row];
   C.delta_shift(adjoint(QShift));

   // solve for the first component
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > ExpIKTol)
   {
      F1[Row] = C;
      LinearSolve(F1[Row], OneMinusTransferRightEA(Op(Row, Row), PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      F1[Row] = C;
      LinearSolve(F1[Row], SubProductRightProjectEA(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}

struct SubProductRightProjectEAOp
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRightProjectEAOp(GenericMPO const& Op_, LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
                              QuantumNumber const& QShift_, MatrixOperator const& Proj_,
                              MatrixOperator const& Ident_, std::complex<double> ExpIK_)
      : Op(Op_), Psi1(Psi1_), Psi2(Psi2_), QShift(QShift_), Proj(Proj_), Ident(Ident_), ExpIK(ExpIK_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& In) const
   {
      MatrixOperator Result = In - ExpIK * delta_shift(inject_right(In, Psi1, Op, Psi2), adjoint(QShift));
      Result -= inner_prod(Proj, Result) * Ident;
      return Result;
   }

   LinearWavefunction const& Psi1;
   LinearWavefunction const& Psi2;
   GenericMPO const& Op;
   QuantumNumber QShift;
   MatrixOperator Proj;
   MatrixOperator Ident;
   std::complex<double> ExpIK;
};

void
SolveStringMPO_EA_Right(MatrixOperator& F1, MatrixOperator const& F0,
                        LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                        LinearWavefunction const& PsiTri,
                        QuantumNumber const& QShift, ProductMPO const& Op,
                        MatrixOperator const& Rho, MatrixOperator const& Ident,
                        std::complex<double> ExpIK, double Tol, int Verbose)
{
   CHECK(Op.is_string());

   if (F1.is_null())
      F1 = MatrixOperator(PsiLeft.Basis2(), PsiRight.Basis2());

   MatrixOperator C = inject_right(F0, PsiTri, Op.data(), PsiRight);
   C.delta_shift(adjoint(QShift));

   if (Rho.is_null() || std::abs(ExpIK - 1.0) > ExpIKTol)
   {
      F1 = C;
      LinearSolve(F1, OneMinusTransferRightEA(Op, PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      F1 = C;
      LinearSolve(F1, SubProductRightProjectEAOp(Op, PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}
