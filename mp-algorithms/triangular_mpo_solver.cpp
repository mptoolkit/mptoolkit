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

#include "triangular_mpo_solver.h"
#include "triangular_mpo_solver_helpers.h"

void
SolveMPO_Left(std::vector<KMatrixPolyType>& EMatK,
              LinearWavefunction const& Psi, QuantumNumber const& QShift,
              BasicTriangularMPO const& Op,
              //std::vector<MatrixOperator> BoundaryE,
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
   EMatK.reserve(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_Left: dimension is " << Dim << std::endl;

   if (EMatK.empty())
   {
      // Make sure the (0,0) part is identity
      DEBUG_TRACE(UnityEpsilon);
      OperatorClassification CheckIdent = classify(Op(0,0), UnityEpsilon);
      if (!CheckIdent.is_identity())
      {
         std::cerr << "SolveMPO_Left: fatal: (0,0) component of the MPO must be the identity operator.\n";
         // the (0,0) component isn't the identity operator, which is a fatal error.
         // Show some diagnosics and quit.
         if (!CheckIdent.is_product())
         {
            std::cerr << "SolveMPO_Left: fatal: component is not a product operator!\n";
         }
         else
         {
            if (CheckIdent.is_unitary())
               std::cerr << "SolveMPO_Left: fatal: component is a unitary operator, but not the identity.\n";
            else if (CheckIdent.is_prop_identity())
            {
               std::cerr << "SolveMPO_Left: fatal: component is proportional "
                         << "to the identity, with prefactor " << CheckIdent.factor() << '\n';
            }
            else
               std::cerr << "SolveMPO_Left: fatal: component has unknown classification.\n";
         }
         PANIC("Fatal")(CheckIdent)(Op(0,0));
      }

      // Initialize the first E matrix.  These are operators acting in the Basis1()
      EMatK.push_back(KMatrixPolyType());
      EMatK[0][1.0] = MatrixPolyType(LeftIdentity);
   }

   int StartCol = EMatK.size();

   // fill out the remainder of EMat with zero
   while (int(EMatK.size()) < Dim)
      EMatK.push_back(KMatrixPolyType());

   // solve recursively column 1 onwards
   for (int Col = StartCol; Col < Dim; ++Col)
   {
      if (Verbose > 0)
      {
         std::cerr << "Solving column " << (Col+1) << " of " << Dim << '\n';
      }

      // Generate the next C matrices, C(n) = sum_{j<Col} Op(j,Col) E_j(n)
      KMatrixPolyType C = inject_left_mask(EMatK, Psi, QShift, Op.data(), Psi, mask_column(Op, Col))[Col];

      // Now do the classification, based on the properties of the diagonal operator
      BasicFiniteMPO Diag = Op(Col, Col);
      OperatorClassification Classification = classify(Diag, UnityEpsilon);

      if (Classification.is_null())
      {
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         if (Verbose > 0)
            std::cerr << "Zero diagonal matrix element at column " << (Col+1) << std::endl;
         EMatK[Col] = SolveZeroDiagonal(C);
      }
      else
      {
         if (Verbose > 0)
            std::cerr << "Non-zero diagonal matrix element at column " << (Col+1) << std::endl;

         // Non-zero diagonal element.
         // In this case we have to consider a possible component with eigenvalue magnitude 1.
         // We call this component (if it exists) the unit matrix.

         KComplexPolyType EParallel;  // components parallel to the identity at some momentum, may be zero

         // Multiplication factor and the left and right eigenvalues.
         std::complex<double> Factor = Classification.factor();
         MatrixOperator UnitMatrixLeft = LeftIdentity;
         MatrixOperator UnitMatrixRight = RightIdentity; //delta_shift(RightIdentity, QShift);

         // If the diagonal operator is unitary, it might have an eigenvalue of magnitude 1
         // Or if we are not in the scalar sector, then there might be an additional eigenvalue 1
         // due to symmetry.  We don't attempt to handle the case where we have more than one
         // eigenvalue 1 in the same symmetry sector.
#if defined(NDEBUG)
         if (Classification.is_unitary() && (!is_scalar(Diag.Basis2()[0]) || !Classification.is_complex_identity()))
#else
         // if we're debugging, then find the eigenvector anyway, even if its the identity
         if (Classification.is_unitary())
#endif
         {
            if (Verbose > 0)
            {
               if (Classification.is_complex_identity())
               {
               std::cerr << "Identity in non-scalar sector at " << Col
                         << ", finding left eigenvalue..." << std::endl;
               }
               else
               {
                  std::cerr << "Non-trivial unitary at column " << Col
                            << ", finding left eigenvalue..." << std::endl;
               }
            }

            // Find the largest eigenvalue
            // We need initial guess vectors in the correct symmetry sector
            UnitMatrixLeft = MakeRandomMatrixOperator(LeftIdentity.Basis1(), LeftIdentity.Basis2(), Diag.Basis2()[0]);

            DEBUG_TRACE(norm_frob(UnitMatrixLeft));
            double lnorm = norm_frob(UnitMatrixLeft);
            //UnitMatrixLeft *= 1.0 / lnorm;
            //UnitMatrixLeft *= 1.0 / norm_frob(UnitMatrixLeft);
            //       UnitMatrixLeft *= 2.0; // adding this brings in spurious components
               std::complex<double> EtaL = FindClosestUnitEigenvalue(UnitMatrixLeft,
                                                                     InjectLeftQShift(Diag, Psi, QShift),
                                                                     Tol, Verbose);
            //UnitMatrixLeft *= lnorm;
            EtaL = std::conj(EtaL); // left eigenvalue, so conjugate (see comment at operator_actions.h)
            if (Verbose > 0)
               std::cerr << "Eigenvalue of unitary operator is " << EtaL << std::endl;
            Factor = EtaL;

            if (std::abs(norm_frob(EtaL) - 1.0) < UnityEpsilon)
            {
               if (Verbose > 0)
                  std::cerr << "Found an eigenvalue 1, so we need the right eigenvector..." << std::endl;

               // we have an eigenvalue of magnitude 1.  Find the right eigenvalue too
               UnitMatrixRight = UnitMatrixLeft;
               UnitMatrixRight = MakeRandomMatrixOperator(LeftIdentity.Basis1(), LeftIdentity.Basis2(), Diag.Basis2()[0]);
               //UnitMatrixRight *= 1.0 / norm_frob(UnitMatrixRight);
               DEBUG_TRACE(norm_frob(UnitMatrixRight));
               double ddd = norm_frob(UnitMatrixRight);
               //UnitMatrixRight *= 1.0 / ddd; //norm_frob(UnitMatrixRight);
               std::complex<double> EtaR = FindClosestUnitEigenvalue(UnitMatrixRight,
                                                                     InjectRightQShift(Diag, Psi,
                                                                                       QShift),
                                                                     Tol, Verbose);
               //UnitMatrixRight *= 3.141;
               if (Verbose > 0)
                  std::cerr << "Right eigenvalue is " << EtaR << std::endl;

               CHECK(norm_frob(EtaL-EtaR) < UnityEpsilon)("Left and right eigenvalues do not agree!")(EtaL)(EtaR);
               // we already determined that the norm is sufficiently close to 1, but
               // fine-tune normalization, which also guarantees that we ultimately set
               // HasEigenvalue1 below
               Factor = EtaL / norm_frob(EtaL);

               //UnitMatrixRight *= 100;

               // normalize the left/right eigenvector pair
               //              UnitMatrixLeft *= 1.0 / (inner_prod(UnitMatrixRight, UnitMatrixLeft));

               //          TRACE(trace(UnitMatrixLeft))(trace(UnitMatrixRight));


               UnitMatrixRight *= 1.0 / (inner_prod(UnitMatrixLeft, UnitMatrixRight));
               CHECK(norm_frob(inner_prod(UnitMatrixLeft, UnitMatrixRight) - 1.0) < 1E-12);
               DEBUG_TRACE(inner_prod(UnitMatrixLeft, UnitMatrixRight));

               //TRACE(trace(UnitMatrixLeft))(trace(UnitMatrixRight));

               //          TRACE(inner_prod(UnitMatrixLeft, LeftIdentity) / norm_frob(UnitMatrixLeft));
               //TRACE(inner_prod(UnitMatrixRight, RightIdentity) / norm_frob(UnitMatrixRight));

            }
            else if (std::abs(norm_frob(EtaL) - 1.0) < UnityEpsilon*100)
            {
               // If we get here, then we have an eigenvalue that is close to 1, but
               // misses the test
               std::cerr << "SolveMPO_Left: warning: Eigenvalue misses tolerance but is near 1"
                         << ", epsilon=" << std::abs(norm_frob(EtaL)-1.0) << '\n';
            }
         }
         else if (Verbose && Classification.is_identity())
         {
            std::cerr << "Diagonal component is the identity\n";
         }
         else if (Verbose && Classification.is_complex_identity())
         {
            std::cerr << "Diagonal component is proportional to the identity\n";
         }

         // If we have an eigenvalue equal to 1, then decompose C into parallel and perpendicular parts
         bool HasEigenvalue1 = false;
         if (std::abs(norm_frob(Factor) - 1.0) < UnityEpsilon)
         {
            HasEigenvalue1 = true;
            //DEBUG_TRACE(UnitMatrixLeft)(UnitMatrixRight);
            if (Verbose > 0)
               std::cerr << "Decomposing parts parallel to the unit matrix\n";
            EParallel = DecomposeParallelPartsWithMomentum(C, Factor, UnitMatrixLeft, UnitMatrixRight, UnityEpsilon, Degree);
         }
         else
         {
            if (Verbose > 0)
            {
               std::cerr << "Diagonal component has spectral radius < 1\n";
               DEBUG_TRACE(Diag)(Classification);
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
            E = DecomposePerpendicularPartsLeft(C, Diag, UnitMatrixLeft, UnitMatrixRight,
                                            Psi, Psi, QShift, 1.0, HasEigenvalue1, Tol, Verbose);
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
               E[I->first][J->first] += std::conj(J->second) * UnitMatrixLeft;
               DEBUG_TRACE(std::conj(J->second));
               DEBUG_TRACE(I->first)(J->first)(inner_prod(E[I->first][J->first], RightIdentity));
            }
         }

         // Finally, set the E matrix element at this column
         //DEBUG_TRACE(E[1.0]);
         EMatK[Col] = E;
      }

   }
}

struct OneMinusTransferLeft2
{
   OneMinusTransferLeft2(GenericMPO const& Op_, LinearWavefunction const& Psi1_,
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

struct SubProductLeftProject2
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeftProject2(LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
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
SolveSimpleMPO_Left2(StateComponent& E1, StateComponent const& E0,
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
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > 1e-14)
   {
      E1[Col] = C;
      LinearSolve(E1[Col], OneMinusTransferLeft2(Op(Col, Col), PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      E1[Col] = C;
      LinearSolve(E1[Col], SubProductLeftProject2(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
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

	 LinearSolve(E1[Col], OneMinusTransferLeft2(Diag, PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
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
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > 1e-14)
   {
      E1[Col] = C;
      LinearSolve(E1[Col], OneMinusTransferLeft2(Op(Col, Col), PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      E1[Col] = C;
      LinearSolve(E1[Col], SubProductLeftProject2(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}

struct SubProductLeftProject2Op
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductLeftProject2Op(GenericMPO const& Op_, LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
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
SolveStringMPO_Left2(MatrixOperator& E1, MatrixOperator const& E0,
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

   if (Rho.is_null() || std::abs(ExpIK - 1.0) > 1e-14)
   {
      E1 = C;
      LinearSolve(E1, OneMinusTransferLeft2(Op, PsiRight, PsiLeft, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      E1 = C;
      LinearSolve(E1, SubProductLeftProject2Op(Op, PsiRight, PsiLeft, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}

struct OneMinusTransferRight2
{
   OneMinusTransferRight2(GenericMPO const& Op_, LinearWavefunction const& Psi1_,
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

struct SubProductRightProject2
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRightProject2(LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
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
SolveSimpleMPO_Right2(StateComponent& F1, StateComponent const& F0,
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
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > 1e-14)
   {
      F1[Row] = C;
      LinearSolve(F1[Row], OneMinusTransferRight2(Op(Row, Row), PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      F1[Row] = C;
      LinearSolve(F1[Row], SubProductRightProject2(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
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

	 LinearSolve(F1[Row], OneMinusTransferRight2(Diag, PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
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
   if (Rho.is_null() || std::abs(ExpIK - 1.0) > 1e-14)
   {
      F1[Row] = C;
      LinearSolve(F1[Row], OneMinusTransferRight2(Op(Row, Row), PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      F1[Row] = C;
      LinearSolve(F1[Row], SubProductRightProject2(PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}

struct SubProductRightProject2Op
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator argument_type;

   SubProductRightProject2Op(GenericMPO const& Op_, LinearWavefunction const& Psi1_, LinearWavefunction const& Psi2_,
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
SolveStringMPO_Right2(MatrixOperator& F1, MatrixOperator const& F0,
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

   if (Rho.is_null() || std::abs(ExpIK - 1.0) > 1e-14)
   {
      F1 = C;
      LinearSolve(F1, OneMinusTransferRight2(Op, PsiLeft, PsiRight, QShift, ExpIK), C, Tol, Verbose);
   }
   else
   {
      // orthogonalize
      C -= inner_prod(Rho, C) * Ident;

      F1 = C;
      LinearSolve(F1, SubProductRightProject2Op(Op, PsiLeft, PsiRight, QShift, Rho, Ident, ExpIK), C, Tol, Verbose);
   }
}
