// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/cg.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "cg.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"

template <typename FwdIter>
double Entropy(FwdIter first, FwdIter last)
{
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue;
      if (x > 0)
         E -= x * log(x);

      ++first;
   }
   return E;
}

template <typename FwdIter>
double TruncError(FwdIter first, FwdIter last)
{
   double E = 0;
   while (first != last)
   {
      double x = first->Eigenvalue;
      if (x > 0)
         E += x;

      ++first;
   }
   return 1.0 - E;
}

Solver::Solver(MPWavefunction const& Psi_, MPOperator const& Op_, MPWavefunction const& Rhs_)
   : x(Psi_), A(Op_), y(Rhs_), Ident(Psi_.GetSymmetryList())
{
   // construct the OpMatrix elements for the right hand side

   typedef MPStateComponent Component;

   // Shift A so the Center matrix lines up with that of x
   while (A.LeftSize() > x.LeftSize()) A.RotateLeft();
   while (A.RightSize() > x.RightSize()) A.RotateRight();

   // Shift y so the Center matrix lines up with that of x
   while (y.LeftSize() > y.LeftSize()) y.RotateLeft();
   while (y.RightSize() > y.RightSize()) y.RotateRight();

   // we could probably do a better job of initializing yprime, if we wanted
   yprime = y;

   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());

   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   yprime_A_x.PushRight(make_vacuum_state(x.GetSymmetryList()));
   yprime_y.PushRight(VacIdent);
   yprime_x.PushRight(VacIdent);

   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      yprime_A_x.PushRight(operator_prod(A.LookupRight(Loc),
                                         yprime.LookupRight(Loc),
                                         yprime_A_x.Right(),
                                         herm(x.LookupRight(Loc))));

      yprime_y.PushRight(operator_prod(yprime.LookupRight(Loc),
                                       yprime_y.Right(),
                                       herm(y.LookupRight(Loc))));

      yprime_x.PushRight(operator_prod(yprime.LookupRight(Loc),
                                       yprime_x.Right(),
                                       herm(x.LookupRight(Loc))));
   }

   //Now start from the left hand side and work inwards
   BasisList LeftVacuum(x.GetSymmetryList());
   CHECK_EQUAL(x.LookupLeft(0).Basis1().size(), 1);
   LeftVacuum.push_back(x.LookupLeft(0).Basis1()[0]);

   VacIdent = MatrixOperator(VectorBasis(LeftVacuum), VectorBasis(LeftVacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   yprime_A_x.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));
   yprime_y.PushLeft(VacIdent);
   yprime_x.PushLeft(VacIdent);

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      yprime_A_x.PushLeft(operator_prod(herm(A.LookupLeft(Loc)),
                                        herm(yprime.LookupLeft(Loc)),
                                        yprime_A_x.Left(),
                                        x.LookupLeft(Loc)));

      yprime_y.PushLeft(operator_prod(herm(yprime.LookupLeft(Loc)),
                                      yprime_y.Left(),
                                      y.LookupLeft(Loc)));

      yprime_x.PushLeft(operator_prod(herm(yprime.LookupLeft(Loc)),
                                      yprime_x.Left(),
                                      x.LookupLeft(Loc)));
   }

   this->DebugCheckBasis();
}

Solver::~Solver()
{
}

void Solver::ExpandLeft()
{
   // correction
   if (!LeftCorrection.is_null())
      LeftCorrection = operator_prod(x.Left(), LeftCorrection, herm(x.Left()));

   // expand x
   x.Left() = prod(x.Left(), x.Center());
   x.Center() = ExpandBasis2(x.Left());
   //TRACE(norm_frob_sq(x.Center()));

   if (!LeftCorrection.is_null())
      LeftCorrection = operator_prod(herm(x.Left()), LeftCorrection, x.Left());

   // for yprime, the Center is y transformed to the basis
     yprime.Left() = prod(yprime.Left(), yprime.Center());
     ExpandBasis2(yprime.Left());

   // matrix elements
   yprime_A_x.PopLeft();
   yprime_A_x.PushLeft(operator_prod(herm(A.Left()),
                                     herm(yprime.Left()),
                                     yprime_A_x.Left(),
                                     x.Left()));

   yprime_y.PopLeft();
   yprime_y.PushLeft(operator_prod(herm(yprime.Left()), yprime_y.Left(), y.Left()));

   yprime_x.PopLeft();
   yprime_x.PushLeft(operator_prod(herm(yprime.Left()), yprime_x.Left(), x.Left()));

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //TRACE(norm_frob_sq(yprime.Center()));

   this->DebugCheckBasis();
}

void Solver::ExpandRight()
{
   //   MatrixOperator Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   //   TRACE(norm_frob(Ax - yprime.Center()));

   // correction
   if (!RightCorrection.is_null())
      RightCorrection = operator_prod(herm(x.Right()), RightCorrection, x.Right());

   // expand x
   x.Right() = prod(x.Center(), x.Right());
   x.Center() = ExpandBasis1(x.Right());

   if (!RightCorrection.is_null())
      RightCorrection = operator_prod(x.Right(), RightCorrection, herm(x.Right()));

#if 0


   yprime_A_x.PopRight();
   yprime_A_x.PushRight(operator_prod(A.Right(),
                                      yprime.Right(),
                                      yprime_A_x.Right(),
                                      herm(x.Right())));

   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   TRACE(norm_frob(Ax - yprime.Center()));





   if (!RightCorrection.is_null())
      RightCorrection = operator_prod(x.Right(), RightCorrection, herm(x.Right()));

   TRACE(yprime.Right().Basis1().total_degree());
   TRACE(yprime.Right().Basis2().total_degree());
   TRACE(trace(scalar_prod(yprime.Right(), herm(yprime.Right()))));
   TRACE(trace(scalar_prod(herm(yprime.Right()), yprime.Right())));

   TRACE(norm_frob(yprime.Center()));
   // for yprime, the Center is y transformed to the basis
   yprime.Right() = prod(yprime.Center(), yprime.Right());

   yprime.Center() = ExpandBasis1(yprime.Right());


   yprime_A_x.PopRight();
   yprime_A_x.PushRight(operator_prod(A.Right(),
                                      yprime.Right(),
                                      yprime_A_x.Right(),
                                      herm(x.Right())));
   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   TRACE(norm_frob(Ax - yprime.Center()));


   DensityMatrix<MatrixOperator> DM(scalar_prod(herm(yprime.Center()), yprime.Center()));
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   while (E != DM.end() && E->Eigenvalue > 1E-14 * DM.EigenSum()) ++E;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   Ax = prod(Ax, adjoint(U), Ax.TransformsAs());
   yprime.Center() = prod(yprime.Center(), adjoint(U), yprime.Center().TransformsAs());
   TRACE(norm_frob(Ax - yprime.Center()));


   yprime.Right() = prod(U, yprime.Right());
   //   yprime.Center() = prod(yprime.Center(), adjoint(U), yprime.Center().TransformsAs());

   //   yprime.Center() = TruncateBasis1(yprime.Right());


   TRACE(norm_frob(yprime.Center()));

   TRACE(yprime.Right().Basis1().total_degree());
   TRACE(yprime.Right().Basis2().total_degree());
   TRACE(trace(scalar_prod(yprime.Right(), herm(yprime.Right()))));
   TRACE(trace(scalar_prod(herm(yprime.Right()), yprime.Right())));

#endif


   yprime.Right() = prod(yprime.Center(), yprime.Right());
   ExpandBasis1(yprime.Right());

   // matrix elements
   yprime_A_x.PopRight();
   yprime_A_x.PushRight(operator_prod(A.Right(),
                                      yprime.Right(),
                                      yprime_A_x.Right(),
                                      herm(x.Right())));



   // yprime_y.PopRight();
   // yprime_y.PushRight(operator_prod(yprime.Right(), yprime_y.Right(), herm(y.Right())));
   //   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   //   TRACE(norm_frob(Ax - yprime.Center()));



   yprime_y.PopRight();
   yprime_y.PushRight(operator_prod(yprime.Right(), yprime_y.Right(), herm(y.Right())));

   yprime_x.PopRight();
   yprime_x.PushRight(operator_prod(yprime.Right(), yprime_x.Right(), herm(x.Right())));

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //TRACE(norm_frob_sq(yprime.Center()));

   this->DebugCheckBasis();
}

void Solver::ShiftRightAndExpand()
{
   MatrixOperator Ax;

   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   TRACE(norm_frob(Ax - yprime.Center()) / norm_frob(yprime.Center()));
   TRACE(norm_frob(Ax));
   TRACE(norm_frob(x));

   // correction
   LeftCorrection = scalar_prod(x.Center(), herm(x.Center()));
   RightCorrection = MatrixOperator();

   // Rotate the wavefunctions to the right.  A and y are fixed and do not
   // need the basis expanded

   //TRACE(norm_frob_sq(yprime.Center()));
   //MatrixOperator Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(),
   // herm(yprime_A_x.Right()));
   //TRACE(norm_frob_sq(Ax));
   //TRACE(norm_frob(Ax - yprime.Center()));

   A.RotateRight();
   y.RotateRight();

   //TRACE(norm_frob_sq(x.Center()));
   x.PushLeft(prod(x.Center(), x.Right()));
   x.PopRight();
   x.Center() = ExpandBasis2(x.Left());
   //TRACE(norm_frob_sq(x.Center()));

   yprime.PushLeft(prod(yprime.Center(), yprime.Right()));
   yprime.PopRight();
   ExpandBasis2(yprime.Left());  // dont update yprime.Center() - do later

   // new matrix elements
   yprime_A_x.PushLeft(operator_prod(herm(A.Left()),
                                     herm(yprime.Left()),
                                     yprime_A_x.Left(),
                                     x.Left()));
   yprime_A_x.PopRight();

   //Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   //TRACE(norm_frob_sq(Ax));

   yprime_y.PushLeft(operator_prod(herm(yprime.Left()), yprime_y.Left(), y.Left()));
   yprime_y.PopRight();

   yprime_x.PushLeft(operator_prod(herm(yprime.Left()), yprime_x.Left(), x.Left()));
   yprime_x.PopRight();

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   TRACE(norm_frob_sq(yprime.Center()));

   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   TRACE(norm_frob(Ax - yprime.Center()) / norm_frob(yprime.Center()));
   TRACE(norm_frob(Ax));
   TRACE(norm_frob(x));

   LeftCorrection = operator_prod(herm(x.Left()), LeftCorrection, x.Left());

   this->DebugCheckBasis();
}

void Solver::ShiftLeftAndExpand()
{
   // correction
   RightCorrection = scalar_prod(herm(x.Center()), x.Center());
   LeftCorrection = MatrixOperator();

   // Rotate the wavefunctions to the left.  A and y are fixed and do not
   // need the basis expanded
   A.RotateLeft();
   y.RotateLeft();

   x.PushRight(prod(x.Left(), x.Center()));
   x.PopLeft();
   x.Center() = ExpandBasis1(x.Right());

   yprime.PushRight(prod(yprime.Left(), yprime.Center()));
   yprime.PopLeft();
   ExpandBasis1(yprime.Right());  // dont update yprime.Center() - do later

   // new matrix elements
   yprime_A_x.PushRight(operator_prod(A.Right(),
                                      yprime.Right(),
                                      yprime_A_x.Right(),
                                      herm(x.Right())));
   yprime_A_x.PopLeft();

   yprime_y.PushRight(operator_prod(yprime.Right(),
                                    yprime_y.Right(),
                                    herm(y.Right())));
   yprime_y.PopLeft();

   yprime_x.PushRight(operator_prod(yprime.Right(),
                                    yprime_x.Right(),
                                    herm(x.Right())));
   yprime_x.PopLeft();

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //TRACE(norm_frob_sq(yprime.Center()));

   RightCorrection = operator_prod(x.Right(), RightCorrection, herm(x.Right()));

   this->DebugCheckBasis();
}

TruncationInfo Solver::TruncateLeft(int MaxStates, double CFactor)
{
   // calculate Ax - TODO: this was already done by the solver, should be saved somewhere
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     yprime_A_x.Left(),
                                     x.Center(),
                                     herm(yprime_A_x.Right()));

   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   TRACE(norm_frob(Ax - yprime.Center()) / norm_frob(yprime.Center()));
   TRACE(norm_frob(Ax));
   TRACE(norm_frob(x));

   MatrixOperator Rho_x = scalar_prod(x.Center(), herm(x.Center()));

   // Calculate the density matrix for yprime
   // This is the sum of the DM for y itself and A*x
   MatrixOperator yp = prod(yprime_y.Left(), y.Center(), y.Center().TransformsAs());
   MatrixOperator Rho_y = scalar_prod(yp, herm(yp));

   MatrixOperator Rho_Ax = operator_prod(yprime_A_x.Left(), Rho_x, herm(yprime_A_x.Left()));

   DensityMatrix<MatrixOperator> RDM((0.5 / trace(Rho_Ax))*Rho_Ax + (0.5 / trace(Rho_y))*Rho_y);

   double const Eps = 1E-4;

#if 0
   if (CFactor != 0)
   {

      TRACE(EigenvaluesHermitian(scalar_prod(herm(yprime_A_x.Left()), yprime_A_x.Left())));

      TRACE(EigenvaluesHermitian(scalar_prod(yprime_A_x.Left(), herm(yprime_A_x.Left()))));

      // Regularize ry = density matrix for y
      MatrixOperator ry = Rho_y;
      MatrixOperator yU = Regularize(ry.Basis1());
      TRACE(trace(Rho_y));
      ry = triple_prod(yU, ry, herm(yU));
      TRACE(trace(ry));

      TRACE(EigenvaluesHermitian(ry));

      ry += Eps * MatrixOperator::make_identity(ry.Basis1());

      TRACE(EigenvaluesHermitian(ry));

      TRACE(trace(ry));
      // Diagonalize ry
      MatrixOperator M = DiagonalizeHermitian(ry);
      TRACE(trace(ry));
      // invert ry
      for (iterator<MatrixOperator>::type I = iterate(ry); I; ++I)
      {
         for (inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
         {
            for (int j = 0; j < ry.Basis1().dim(J.index1()); ++j)
            {
               if (LinearAlgebra::norm_2((*J)(j,j)) < 1E-16)
                  (*J)(j,j) = 0;
               else
                  (*J)(j,j) = 1.0 / (*J)(j,j);
            }
         }
      }
      TRACE(trace(ry));
      TRACE(EigenvaluesHermitian(ry));

       // convert back to our basis
      ry = triple_prod(herm(M), ry, M);
      TRACE(trace(ry));
      ry = triple_prod(herm(yU), ry, yU);
      TRACE(trace(ry));
      // apply the A matrices
      MatrixOperator Rho_ry = operator_prod(herm(yprime_A_x.Left()), ry, yprime_A_x.Left());
      // regularize again so we can invert
      MatrixOperator iU = Regularize(Rho_ry.Basis1());
      Rho_ry = triple_prod(iU, Rho_ry, herm(iU));
      DEBUG_TRACE(trace(Rho_ry));

      TRACE(EigenvaluesHermitian(Rho_ry));

      // Diagonalize and invert
      M = DiagonalizeHermitian(Rho_ry);
      DEBUG_TRACE(trace(Rho_ry));
      //      Rho_ry += Eps * MatrixOperator::make_identity(Rho_ry.Basis1());
      for (iterator<MatrixOperator>::type I = iterate(Rho_ry); I; ++I)
      {
         for (inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
         {
            for (int j = 0; j < Rho_ry.Basis1().dim(J.index1()); ++j)
            {
               if (LinearAlgebra::norm_2((*J)(j,j)) < 1E-5)
                  (*J)(j,j) = 0;
               else
                  (*J)(j,j) = 1.0 / (*J)(j,j);
            }
         }
      }

      TRACE(EigenvaluesHermitian(Rho_ry));

      DEBUG_TRACE(trace(Rho_ry));
      // convert back to our basis
      Rho_ry = triple_prod(herm(M), Rho_ry, M);
      Rho_ry = triple_prod(herm(iU), Rho_ry, iU);

      // Add the correction
      DEBUG_TRACE(trace(Rho_ry));
      Rho_x += (CFactor / trace(Rho_ry)) * Rho_ry;
   }
#else
   if (!LeftCorrection.is_null() && CFactor != 0)
      Rho_x += CFactor * LeftCorrection;
#endif

   // Calculate the density matrix for x

   MatrixOperator Rho_xy = operator_prod(herm(yprime_A_x.Left()), Rho_y, yprime_A_x.Left());

   DensityMatrix<MatrixOperator> DM((0.5 / trace(Rho_x))*Rho_x);// + (0.5/trace(Rho_xy))*Rho_xy);
   //   DM.DensityMatrixReport(std::cerr);
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
      //          (E->Eigenvalue > EigenvalueEpsilon * DM.EigenSum()
      //           || Count < MaxStates))
   {
      ++E; ++Count;
   }

   TruncationInfo Result;
   Result.m = Count;
   Result.entropy = Entropy(DM.begin(), E);
   Result.trunc = TruncError(DM.begin(), E);




   // truncate
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   x.Center() = prod(U, x.Center(), Ident);
   x.Left() = prod(x.Left(), adjoint(U));

   yprime_A_x.Left() = prod(yprime_A_x.Left(), adjoint(U));
   yprime_x.Left() = prod(yprime_x.Left(), adjoint(U), Ident);

   //DensityMatrix<MatrixOperator> RDM(0.5*Rho_Ax + 0.5*Rho_y);

   //   DM.DensityMatrixReport(std::cerr);
   E = RDM.begin();
   Count = 0;
   while (E != RDM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   //   //TRACE(Count);
   U = RDM.ConstructTruncator(RDM.begin(), E);

   // truncate
   yprime.Center() = prod(U, yprime.Center(), Ident);
   yprime.Left() = prod(yprime.Left(), adjoint(U));

   yprime_A_x.Left() = prod(U, yprime_A_x.Left());
   yprime_y.Left() = prod(U, yprime_y.Left(), Ident);
   yprime_x.Left() = prod(U, yprime_x.Left(), Ident);


   Ax = operator_prod(A.Center(), yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
   TRACE(norm_frob(Ax - yprime.Center()) / norm_frob(yprime.Center()));
   TRACE(norm_frob(Ax));
   TRACE(norm_frob(x));



   this->DebugCheckBasis();
   return Result;
}

TruncationInfo Solver::TruncateRight(int MaxStates, double CFactor)
{
   // calculate Ax - TODO: this was already done by the solver, should be saved somewhere
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     yprime_A_x.Left(),
                                     x.Center(),
                                     herm(yprime_A_x.Right()));

   MatrixOperator Rho_x = scalar_prod(herm(x.Center()), x.Center());

   // Calculate the density matrix for yprime
   // This is the sum of the DM for y itself and A*x
   // It is perhaps better to use the truncated x to calculate Ax here
   MatrixOperator Rho_Ax = scalar_prod(herm(Ax), Ax);
   MatrixOperator Rho_y = scalar_prod(herm(yprime.Center()), yprime.Center());

   MatrixOperator yp = prod(y.Center(), adjoint(yprime_y.Right()), y.Center().TransformsAs());
   Rho_y = scalar_prod(herm(yp), yp);

   Rho_Ax = operator_prod(yprime_A_x.Right(), Rho_x, herm(yprime_A_x.Right()));

   DensityMatrix<MatrixOperator> RDM((0.5 / trace(Rho_Ax))*Rho_Ax + (0.5 / trace(Rho_y))*Rho_y);



   //   if (!RightCorrection.is_null() && CFactor != 0)
   //      Rho_x += CFactor * RightCorrection;

   // Calculate the density matrix for x
   MatrixOperator Rho_xy = operator_prod(herm(yprime_A_x.Right()), Rho_y, yprime_A_x.Right());

   DensityMatrix<MatrixOperator> DM((0.5 / trace(Rho_x))*Rho_x);// + (0.5/trace(Rho_xy))*Rho_xy);
   //   DM.DensityMatrixReport(std::cerr);
   DensityMatrix<MatrixOperator>::const_iterator E = DM.begin();
   int Count = 0;
   while (E != DM.end() && Count < MaxStates)
      //          (E->Eigenvalue > EigenvalueEpsilon * DM.EigenSum()
      //           || Count < MaxStates))
   {
      ++E; ++Count;
   }

   TruncationInfo Result;
   Result.m = Count;
   Result.entropy = Entropy(DM.begin(), E);
   Result.trunc = TruncError(DM.begin(), E);





   // truncate
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   x.Center() = prod(x.Center(), adjoint(U), Ident);
   x.Right() = prod(U, x.Right());

   yprime_A_x.Right() = prod(yprime_A_x.Right(), adjoint(U));
   yprime_x.Right() = prod(yprime_x.Right(), adjoint(U), Ident);

   E = RDM.begin();
   Count = 0;
   while (E != RDM.end() && Count < MaxStates)
   {
      ++E;
      ++Count;
   }
   //TRACE(Count);
   U = RDM.ConstructTruncator(RDM.begin(), E);
   // truncate

   //   //TRACE(prod(U, adjoint(U), Ident));

   yprime.Center() = prod(yprime.Center(), adjoint(U), Ident);
   yprime.Right() = prod(U, yprime.Right());

   yprime_A_x.Right() = prod(U, yprime_A_x.Right());
   yprime_y.Right() = prod(U, yprime_y.Right(), Ident);
   yprime_x.Right() = prod(U, yprime_x.Right(), Ident);

   this->DebugCheckBasis();
   return Result;
}

void Solver::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(yprime_A_x.Left().Basis1(), yprime.Center().Basis1());
   DEBUG_CHECK_EQUAL(yprime_A_x.Right().Basis1(), yprime.Center().Basis2());

   DEBUG_CHECK_EQUAL(yprime_A_x.Left().Basis2(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(yprime_A_x.Right().Basis2(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(yprime_y.Left().Basis1(), yprime.Center().Basis1());
   DEBUG_CHECK_EQUAL(yprime_y.Right().Basis1(), yprime.Center().Basis2());

   DEBUG_CHECK_EQUAL(yprime_y.Left().Basis2(), y.Center().Basis1());
   DEBUG_CHECK_EQUAL(yprime_y.Right().Basis2(), y.Center().Basis2());

   DEBUG_CHECK_EQUAL(yprime_x.Left().Basis1(), yprime.Center().Basis1());
   DEBUG_CHECK_EQUAL(yprime_x.Right().Basis1(), yprime.Center().Basis2());

   DEBUG_CHECK_EQUAL(yprime_x.Left().Basis2(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(yprime_x.Right().Basis2(), x.Center().Basis2());
}
