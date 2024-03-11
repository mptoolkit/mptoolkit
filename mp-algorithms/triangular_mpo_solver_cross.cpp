// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/triangular_mpo_solver_cross.cpp
//
// Copyright (C) 2009-2022 Ian McCulloch <ian@qusim.net>
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
SolveMPO_Left_Cross(std::vector<KMatrixPolyType>& EMatK,
                    LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                    BasicTriangularMPO const& Op, MatrixOperator const& LeftIdentity,
                    MatrixOperator const& RightIdentity, std::complex<double> lambda, bool NeedFinalMatrix,
                    int Degree, double Tol,
                    double UnityEpsilon, int Verbose)
{
   CHECK_EQUAL(RightIdentity.Basis1(), Psi1.Basis1());
   CHECK_EQUAL(RightIdentity.Basis2(), Psi2.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis1(), Psi1.Basis1());
   CHECK_EQUAL(LeftIdentity.Basis2(), Psi2.Basis1());

   DEBUG_TRACE(Verbose)(Degree)(Tol);

   // lambda is the leading eigenvalue of the transfer matrix.  We solve the MPO with respect to
   // <Psi1|Op|Psi2> / <Psi1|Psi2>.  This amounts to dividing through by lambda every time we contract over a unit cell.
   // Or equivalently, scale everything by Scale ( = 1 / lambda).
   std::complex<double> Scale = 1.0 / lambda;

   int Dim = Op.Basis1().size();       // dimension of the MPO
   EMatK.reserve(Dim);

   if (Verbose > 0)
      std::cerr << "SolveMPO_Left_Cross: dimension is " << Dim << std::endl;

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
      KMatrixPolyType C = inject_left_mask(EMatK, Psi1, QShift, Op.data(), Psi2, mask_column(Op, Col))[Col];
      ScalePoly(C, Scale);

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
                                                                  InjectLeftQShift(Psi1, QShift, Diag, Psi2, Scale),
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
                                                                     InjectRightQShift(Psi1, QShift, Diag, Psi2, Scale),
                                                                     Tol, Verbose);
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
                                            Psi1, Psi2, QShift, Scale, HasEigenvalue1, Tol, Verbose);
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
