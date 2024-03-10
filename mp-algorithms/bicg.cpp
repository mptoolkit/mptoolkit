// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/bicg.cpp
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

/*
  2005-10-12: NOTE: this header is now obsolete.
*/

#warning Obsolete

#include "bicg.h"
#include "biconjugategradient.h"

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

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiply(SimpleOperator const& Op_,
                      MPStateComponent const& Left_,
                      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, Left, Psi, herm(Right));
   }

   SimpleOperator const& Op;
   MPStateComponent const& Left;
   MPStateComponent const& Right;
};

inline
SuperblockMultiply::SuperblockMultiply(SimpleOperator const& Op_,
                                       MPStateComponent const& Left_,
                                       MPStateComponent const& Right_)
   : Op(Op_), Left(Left_), Right(Right_)
{
}

struct SuperblockMultiplyHerm
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiplyHerm(SimpleOperator const& Op_,
                          MPStateComponent const& Left_,
                          MPStateComponent const& Right_)
      : Op(Op_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, herm(Left), Psi, Right);
   }

   SimpleOperator const& Op;
   MPStateComponent const& Left;
   MPStateComponent const& Right;
};

struct MapRhsToLhs
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   MapRhsToLhs(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      return triple_prod(herm(Left), r, Right);
   }

   MatrixOperator const& Left;
   MatrixOperator const& Right;
};

struct MapLhsToRhs
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   MapLhsToRhs(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      return triple_prod(Left, r, herm(Right));
   }

   MatrixOperator const& Left;
   MatrixOperator const& Right;
};

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

double Solver::Energy()
{
   //   DEBUG_TRACE(OpMatrices.Left())(OpMatrices.Right());
   //   DEBUG_TRACE(scalar_prod(Psi.Center(), Psi.Center()));
   //   DEBUG_TRACE(scalar_prod(PsiP, PsiP));

   return inner_prod(operator_prod(conj(A.Center()),
                                   yprime_A_x.Left(),
                                   x.Center(),
                                   herm(yprime_A_x.Right())),
                     yprime.Center()).real();
}

double Solver::Solve(int MaxIterations)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   double Tol = 1E-10;

   BiConjugateGradient(x.Center(),
                       SuperblockMultiply(conj(A.Center()),
                                          yprime_A_x.Left(),
                                          yprime_A_x.Right()),
                       SuperblockMultiplyHerm(conj(A.Center()),
                                              yprime_A_x.Left(),
                                              yprime_A_x.Right()),
                       yprime.Center(),
                       Iter, Tol,
                       MapRhsToLhs(yprime_x.Left(), yprime_x.Right()),
                       MapLhsToRhs(yprime_x.Left(), yprime_x.Right()));


   //   TRACE(Iter)(Tol);
   return Tol;
}

void Solver::Evaluate()
{
   yprime.Center() = operator_prod(conj(A.Center()),
                                   yprime_A_x.Left(), x.Center(), herm(yprime_A_x.Right()));
}

void Solver::ExpandLeft()
{
   // expand x
   x.Left() = prod(x.Left(), x.Center());
   x.Center() = ExpandBasis2(x.Left());

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
   yprime_y.PushLeft(operator_prod(herm(yprime.Left()),
                                   yprime_y.Left(),
                                   y.Left()));

   yprime_x.PopLeft();
   yprime_x.PushLeft(operator_prod(herm(yprime.Left()),
                                   yprime_x.Left(),
                                   x.Left()));

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //   TRACE(norm_frob_sq(yprime.Center()));

   this->DebugCheckBasis();
}

void Solver::ExpandRight()
{
   // expand x
   x.Right() = prod(x.Center(), x.Right());
   x.Center() = ExpandBasis1(x.Right());

   // for yprime, the Center is y transformed to the basis
   yprime.Right() = prod(yprime.Center(), yprime.Right());
   ExpandBasis1(yprime.Right());

   // matrix elements
   yprime_A_x.PopRight();
   yprime_A_x.PushRight(operator_prod(A.Right(),
                                      yprime.Right(),
                                      yprime_A_x.Right(),
                                      herm(x.Right())));

   yprime_y.PopRight();
   yprime_y.PushRight(operator_prod(yprime.Right(), yprime_y.Right(), herm(y.Right())));

   yprime_x.PopRight();
   yprime_x.PushRight(operator_prod(yprime.Right(), yprime_x.Right(), herm(x.Right())));

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //   TRACE(norm_frob_sq(yprime.Center()));

   this->DebugCheckBasis();
}

void Solver::ShiftRightAndExpand()
{
   // Rotate the wavefunctions to the right.  A and x are fixed and do not
   // need the basis expanded
   A.RotateRight();
   y.RotateRight();

   x.PushLeft(prod(x.Center(), x.Right()));
   x.PopRight();
   x.Center() = ExpandBasis2(x.Left());

   yprime.PushLeft(prod(yprime.Center(), yprime.Right()));
   yprime.PopRight();
   ExpandBasis2(yprime.Left());  // dont update yprime.Center() - do later

   // new matrix elements
   yprime_A_x.PushLeft(operator_prod(herm(A.Left()), herm(yprime.Left()), yprime_A_x.Left(), x.Left()));
   yprime_A_x.PopRight();

   yprime_y.PushLeft(operator_prod(herm(yprime.Left()), yprime_y.Left(), y.Left()));
   yprime_y.PopRight();

   yprime_x.PushLeft(operator_prod(herm(yprime.Left()), yprime_x.Left(), x.Left()));
   yprime_x.PopRight();

   // yprime is just y transformed into the new basis
   yprime.Center() = triple_prod(yprime_y.Left(), y.Center(), herm(yprime_y.Right()));
   //   TRACE(norm_frob_sq(yprime.Center()));

   this->DebugCheckBasis();
}

void Solver::ShiftLeftAndExpand()
{
   // Rotate the wavefunctions to the left.  A and x are fixed and do not
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
   //   TRACE(norm_frob_sq(yprime.Center()));

   this->DebugCheckBasis();
}

TruncationInfo Solver::TruncateLeft(int MaxStates)
{
   // calculate Ax - TODO: this was already done by the solver, should be saved somewhere
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     yprime_A_x.Left(),
                                     x.Center(),
                                     herm(yprime_A_x.Right()));

    MatrixOperator Rho_x = scalar_prod(x.Center(), herm(x.Center()));
   // Calculate the density matrix for yprime
   // This is the sum of the DM for y itself and A*x
   MatrixOperator Rho_Ax = scalar_prod(Ax, herm(Ax));
   MatrixOperator Rho_y = scalar_prod(yprime.Center(), herm(yprime.Center()));

   MatrixOperator yp = prod(yprime_y.Left(), y.Center(), y.Center().TransformsAs());
   Rho_y = scalar_prod(yp, herm(yp));

   Rho_Ax = operator_prod(yprime_A_x.Left(), Rho_x, herm(yprime_A_x.Left()));


  // Calculate the density matrix for x
   DensityMatrix<MatrixOperator> DM(Rho_x);

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

   DensityMatrix<MatrixOperator> RDM((0.5 / trace(Rho_Ax))*Rho_Ax + (0.5 / trace(Rho_y))*Rho_y);

   //   DM.DensityMatrixReport(std::cerr);
   E = RDM.begin();
   Count = 0;
   while (E != RDM.end() && Count < MaxStates)
   {
      ++E; ++Count;
   }
   //   TRACE(Count);
   U = RDM.ConstructTruncator(RDM.begin(), E);

   // truncate
   yprime.Center() = prod(U, yprime.Center(), Ident);
   yprime.Left() = prod(yprime.Left(), adjoint(U));

   yprime_A_x.Left() = prod(U, yprime_A_x.Left());
   yprime_y.Left() = prod(U, yprime_y.Left(), Ident);
   yprime_x.Left() = prod(U, yprime_x.Left(), Ident);

   this->DebugCheckBasis();
   return Result;
}

TruncationInfo Solver::TruncateRight(int MaxStates)
{
   // calculate Ax - TODO: this was already done by the solver, should be saved somewhere
   MatrixOperator Ax = operator_prod(conj(A.Center()),
                                     yprime_A_x.Left(),
                                     x.Center(),
                                     herm(yprime_A_x.Right()));

   MatrixOperator Rho_x = scalar_prod(herm(x.Center()), x.Center());

   MatrixOperator Rho_Ax = scalar_prod(herm(Ax), Ax);
   MatrixOperator Rho_y = scalar_prod(herm(yprime.Center()), yprime.Center());

   MatrixOperator yp = prod(y.Center(), adjoint(yprime_y.Right()), y.Center().TransformsAs());
   Rho_y = scalar_prod(herm(yp), yp);

   Rho_Ax = operator_prod(yprime_A_x.Right(), Rho_x, herm(yprime_A_x.Right()));

   // Calculate the density matrix for x
   DensityMatrix<MatrixOperator> DM(Rho_x);
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

   // Calculate the density matrix for yprime
   // This is the sum of the DM for y itself and A*x
   DensityMatrix<MatrixOperator> RDM((0.5 / trace(Rho_Ax))*Rho_Ax + (0.5 / trace(Rho_y))*Rho_y);

   E = RDM.begin();
   Count = 0;
   while (E != RDM.end() && Count < MaxStates)
   {
      ++E;
      ++Count;
   }
   TRACE(Count);
   U = RDM.ConstructTruncator(RDM.begin(), E);
   // truncate

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
