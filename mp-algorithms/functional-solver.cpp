// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/functional-solver.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "functional-solver.h"
#include "ddmrg_functional.h"
#include "gmres.h"

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiply(SimpleOperator const& Op_,
                      MPStateComponent const& Left_,
                      MPStateComponent const& Right_)
      : Op(Op_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, Left, Psi, herm(Right));
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
};

struct SuperblockMultiplySquare
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiplySquare(SimpleOperator const& Op_,
                            MPStateComponent const& Left_,
                            MPStateComponent const& Right_,
                            SimpleOperator const& Op2_,
                            MPStateComponent const& Left2_,
                            MPStateComponent const& Right2_)
      : Op(Op_), Left(Left_), Right(Right_), Op2(Op2_), Left2(Left2_), Right2(Right2_){}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op2, Left2, operator_prod(Op, Left, Psi, herm(Right)), herm(Right2));
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
   SimpleOperator Op2;
   MPStateComponent Left2, Right2;
};

FunctionalSolver::FunctionalSolver(CenterWavefunction const& Psi_, SplitOperator const& A_,
               SplitOperator const& B_, SplitOperator const& H_, SplitOperator const& Hb_,
               CenterWavefunction const& Rhs_, double Freq_, double Broad_)
   : x(Psi_), A(A_), B(B_), H(H_), Hb(Hb_), y(Rhs_), Frequency(Freq_),
     Broadening(Broad_), Ident(Psi_.GetSymmetryList()), Precision(0.0), MinIterations(0), TwoStepKrylov(false),
     MixNormalize(false), LanczosMixFactor(0.0), LinearSolver(false), SquareMeHarder(false), UseResid(false)
{
   // construct the OpMatrix elements for the right hand side

   typedef MPStateComponent Component;

   // Shift A so the Center matrix lines up with that of x
   while (A.LeftSize() > x.LeftSize()) A.RotateLeft();
   while (A.RightSize() > x.RightSize()) A.RotateRight();

   while (B.LeftSize() > x.LeftSize()) B.RotateLeft();
   while (B.RightSize() > x.RightSize()) B.RotateRight();

   while (H.LeftSize() > x.LeftSize()) H.RotateLeft();
   while (H.RightSize() > x.RightSize()) H.RotateRight();

   while (Hb.LeftSize() > x.LeftSize()) Hb.RotateLeft();
   while (Hb.RightSize() > x.RightSize()) Hb.RotateRight();

   // Shift y so the Center matrix lines up with that of x
   while (y.LeftSize() > y.LeftSize()) y.RotateLeft();
   while (y.RightSize() > y.RightSize()) y.RotateRight();


   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());

   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   x_A_x.PushRight(make_vacuum_state(x.GetSymmetryList()));
   x_H_x.PushRight(make_vacuum_state(x.GetSymmetryList()));
   x_Hb_x.PushRight(make_vacuum_state(x.GetSymmetryList()));
   x_B_y.PushRight(make_vacuum_state(x.GetSymmetryList()));
   x_y.PushRight(VacIdent);

   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      x_A_x.PushRight(operator_prod(A.LookupRight(Loc),
                                    x.LookupRight(Loc),
                                    x_A_x.Right(),
                                    herm(x.LookupRight(Loc))));

      x_H_x.PushRight(operator_prod(H.LookupRight(Loc),
                                    x.LookupRight(Loc),
                                    x_H_x.Right(),
                                    herm(x.LookupRight(Loc))));

      x_Hb_x.PushRight(operator_prod(Hb.LookupRight(Loc),
                                     x.LookupRight(Loc),
                                     x_Hb_x.Right(),
                                     herm(x.LookupRight(Loc))));

      x_B_y.PushRight(operator_prod(B.LookupRight(Loc),
                                    x.LookupRight(Loc),
                                    x_B_y.Right(),
                                    herm(y.LookupRight(Loc))));

      x_y.PushRight(operator_prod(x.LookupRight(Loc),
                                  x_y.Right(),
                                  herm(y.LookupRight(Loc))));
   }

   //Now start from the left hand side and work inwards
   BasisList LeftVacuum(x.GetSymmetryList());
   CHECK_EQUAL(x.LookupLeft(0).Basis1().size(), 1);
   LeftVacuum.push_back(x.LookupLeft(0).Basis1()[0]);

   VacIdent = MatrixOperator(VectorBasis(LeftVacuum), VectorBasis(LeftVacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   x_A_x.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));
   x_H_x.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));
   x_Hb_x.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));
   x_B_y.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));
   x_y.PushLeft(VacIdent);

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      x_A_x.PushLeft(operator_prod(herm(A.LookupLeft(Loc)),
                                   herm(x.LookupLeft(Loc)),
                                   x_A_x.Left(),
                                   x.LookupLeft(Loc)));

      x_H_x.PushLeft(operator_prod(herm(H.LookupLeft(Loc)),
                                   herm(x.LookupLeft(Loc)),
                                   x_H_x.Left(),
                                   x.LookupLeft(Loc)));

      x_Hb_x.PushLeft(operator_prod(herm(Hb.LookupLeft(Loc)),
                                    herm(x.LookupLeft(Loc)),
                                    x_Hb_x.Left(),
                                    x.LookupLeft(Loc)));

      x_B_y.PushLeft(operator_prod(herm(B.LookupLeft(Loc)),
                                   herm(x.LookupLeft(Loc)),
                                   x_B_y.Left(),
                                   y.LookupLeft(Loc)));

      x_y.PushLeft(operator_prod(herm(x.LookupLeft(Loc)),
                                 x_y.Left(),
                                 y.LookupLeft(Loc)));
   }

   this->DebugCheckBasis();
}

FunctionalSolver::~FunctionalSolver()
{
}

void FunctionalSolver::ExpandLeft()
{
   // expand x
   x.Left() = prod(x.Left(), x.Center());
   x.Center() = ExpandBasis2(x.Left());

   // matrix elements
   x_A_x.PopLeft();
   x_A_x.PushLeft(operator_prod(herm(A.Left()),
                                herm(x.Left()),
                                x_A_x.Left(),
                                x.Left()));

   x_H_x.PopLeft();
   x_H_x.PushLeft(operator_prod(herm(H.Left()),
                                herm(x.Left()),
                                x_H_x.Left(),
                                x.Left()));

   x_Hb_x.PopLeft();
   x_Hb_x.PushLeft(operator_prod(herm(Hb.Left()),
                                 herm(x.Left()),
                                 x_Hb_x.Left(),
                                 x.Left()));

   x_B_y.PopLeft();
   x_B_y.PushLeft(operator_prod(herm(B.Left()),
                                herm(x.Left()),
                                x_B_y.Left(),
                                y.Left()));

   x_y.PopLeft();
   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));

   this->DebugCheckBasis();
}

void FunctionalSolver::ExpandRight()
{
   // expand x
   x.Right() = prod(x.Center(), x.Right());
   x.Center() = ExpandBasis1(x.Right());

   // matrix elements
   x_A_x.PopRight();
   x_A_x.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_x.Right(),
                                 herm(x.Right())));

   x_H_x.PopRight();
   x_H_x.PushRight(operator_prod(H.Right(),
                                 x.Right(),
                                 x_H_x.Right(),
                                 herm(x.Right())));

   x_Hb_x.PopRight();
   x_Hb_x.PushRight(operator_prod(Hb.Right(),
                                  x.Right(),
                                  x_Hb_x.Right(),
                                  herm(x.Right())));

   x_B_y.PopRight();
   x_B_y.PushRight(operator_prod(B.Right(),
                                 x.Right(),
                                 x_B_y.Right(),
                                 herm(y.Right())));

   x_y.PopRight();
   x_y.PushRight(operator_prod(x.Right(), x_y.Right(), herm(y.Right())));

   this->DebugCheckBasis();
}

void FunctionalSolver::ShiftRightAndExpand()
{
   // Rotate the wavefunctions to the right.  A and y are fixed and do not
   // need the basis expanded

   A.RotateRight();
   H.RotateRight();
   Hb.RotateRight();
   B.RotateRight();
   y.RotateRight();

   //TRACE(norm_frob_sq(x.Center()));
   x.PushLeft(prod(x.Center(), x.Right()));
   x.PopRight();
   x.Center() = ExpandBasis2(x.Left());
   //TRACE(norm_frob_sq(x.Center()));

   // new matrix elements
   x_A_x.PushLeft(operator_prod(herm(A.Left()),
                                herm(x.Left()),
                                x_A_x.Left(),
                                x.Left()));
   x_A_x.PopRight();

   x_H_x.PushLeft(operator_prod(herm(H.Left()),
                                herm(x.Left()),
                                x_H_x.Left(),
                                x.Left()));
   x_H_x.PopRight();

   x_Hb_x.PushLeft(operator_prod(herm(Hb.Left()),
                                 herm(x.Left()),
                                 x_Hb_x.Left(),
                                 x.Left()));
   x_Hb_x.PopRight();

   x_B_y.PushLeft(operator_prod(herm(B.Left()),
                                herm(x.Left()),
                                x_B_y.Left(),
                                y.Left()));
   x_B_y.PopRight();

   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));
   x_y.PopRight();

   this->DebugCheckBasis();
}

void FunctionalSolver::ShiftLeftAndExpand()
{
   // Rotate the wavefunctions to the left.  A and y are fixed and do not
   // need the basis expanded
   A.RotateLeft();
   H.RotateLeft();
   Hb.RotateLeft();
   B.RotateLeft();
   y.RotateLeft();

   x.PushRight(prod(x.Left(), x.Center()));
   x.PopLeft();
   x.Center() = ExpandBasis1(x.Right());

   // new matrix elements
   x_A_x.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_x.Right(),
                                 herm(x.Right())));
   x_A_x.PopLeft();

   x_H_x.PushRight(operator_prod(H.Right(),
                                 x.Right(),
                                 x_H_x.Right(),
                                 herm(x.Right())));
   x_H_x.PopLeft();

   x_Hb_x.PushRight(operator_prod(Hb.Right(),
                                  x.Right(),
                                  x_Hb_x.Right(),
                                  herm(x.Right())));
   x_Hb_x.PopLeft();

   x_B_y.PushRight(operator_prod(B.Right(),
                                 x.Right(),
                                 x_B_y.Right(),
                                 herm(y.Right())));
   x_B_y.PopLeft();

   x_y.PushRight(operator_prod(x.Right(),
                               x_y.Right(),
                               herm(y.Right())));
   x_y.PopLeft();

   this->DebugCheckBasis();
}

TruncationInfo FunctionalSolver::TruncateLeft(StatesInfo const& SInfo, double CFactor)
{
   // Calculate the density matrix for x
   MatrixOperator yp = prod(x_y.Left(), y.Center(), y.Center().TransformsAs());
   MatrixOperator Rho_y = scalar_prod(yp, herm(yp));
   MatrixOperator Rho_x = scalar_prod(x.Center(), herm(x.Center()));

   if (MixNormalize)
   {
      Rho_x *= 1.0 / trace(Rho_x);
      Rho_y *= 1.0 / trace(Rho_y);
   }

   MatrixOperator Rho = (1.0 / (1.0 + LanczosMixFactor)) * Rho_x
      + (LanczosMixFactor / (1.0 + LanczosMixFactor)) * Rho_y;

   if (CFactor != 0)
   {
      for (unsigned i = 0; i < x_A_x.Left().size(); ++i)
      {
         MatrixOperator Correction = triple_prod(x_A_x.Left()[i],
                                                 Rho,
                                                 herm(x_A_x.Left()[i]));
         Correction *= (CFactor / trace(Correction) * double(x_A_x.Left().size()));
         Rho += Correction;
      }
#if 0
      MatrixOperator Rho_y = scalar_prod(y.Center(), herm(y.Center()));
      Rho_y = triple_prod(x_y.Left(), Rho_y, herm(x_y.Left()));
      MatrixOperator Correction = operator_prod(x_A_x.Left(), Rho_y, herm(x_A_x.Left()));
      MatrixOperator xAxInv = scalar_prod(x_A_x.Left(), herm(x_A_x.Left()));
      Tensor::InvertIrregularHPD(xAxInv);
      Correction = triple_prod(xAxInv, Correction, herm(xAxInv));
      DensityMat += (CFactor / trace(Correction)) * Correction;
#endif
   }

   DensityMatrix<MatrixOperator> DM(Rho);
   //   DM.DensityMatrixReport(std::cerr);
   TruncationInfo Info;
   DensityMatrix<MatrixOperator>::const_iterator E
      = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), SInfo, Info);

   // truncate
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   x.Center() = prod(U, x.Center(), Ident);
   x.Left() = prod(x.Left(), adjoint(U));

   x_A_x.Left() = prod(prod(U, x_A_x.Left()), adjoint(U));
   x_H_x.Left() = prod(prod(U, x_H_x.Left()), adjoint(U));
   x_Hb_x.Left() = prod(prod(U, x_Hb_x.Left()), adjoint(U));
   x_B_y.Left() = prod(U, x_B_y.Left());
   x_y.Left() = prod(U, x_y.Left(), Ident);

   this->DebugCheckBasis();
   return Info;
}

TruncationInfo FunctionalSolver::TruncateRight(StatesInfo const& SInfo, double CFactor)
{
   MatrixOperator yp = y.Center() * herm(x_y.Right());
   MatrixOperator Rho_y = scalar_prod(herm(yp), yp);
   MatrixOperator Rho_x = scalar_prod(herm(x.Center()), x.Center());

   if (MixNormalize)
   {
      Rho_x *= 1.0 / trace(Rho_x);
      Rho_y *= 1.0 / trace(Rho_y);
   }

   MatrixOperator Rho = (1.0 / (1.0 + LanczosMixFactor)) * Rho_x
      + (LanczosMixFactor / (1.0 + LanczosMixFactor)) * Rho_y;

   if (CFactor != 0)
   {
      for (unsigned i = 0; i < x_A_x.Right().size(); ++i)
      {
         MatrixOperator Correction = triple_prod(x_A_x.Right()[i],
                                                 Rho,
                                                 herm(x_A_x.Right()[i]));
         Correction *= CFactor / trace(Correction) * double(x_A_x.Right().size());
         Rho += Correction;
      }
#if 0
      MatrixOperator Rho_y = scalar_prod(herm(y.Center()), y.Center());
      Rho_y = triple_prod(x_y.Right(), Rho_y, herm(x_y.Right()));
      MatrixOperator Correction = operator_prod(x_A_x.Right(), Rho_y, herm(x_A_x.Right()));
      MatrixOperator xAxInv = scalar_prod(x_A_x.Right(), herm(x_A_x.Right()));
      Tensor::InvertIrregularHPD(xAxInv);
      Correction = triple_prod(xAxInv, Correction, herm(xAxInv));
      DensityMat += (CFactor / trace(Correction)) * Correction;
#endif
   }

   DensityMatrix<MatrixOperator> DM(Rho);
   //   DM.DensityMatrixReport(std::cerr);
   TruncationInfo Info;
   DensityMatrix<MatrixOperator>::const_iterator E
      = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), SInfo, Info);

   // truncate
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   x.Center() = prod(x.Center(), adjoint(U), Ident);
   x.Right() = prod(U, x.Right());

   x_A_x.Right() = prod(prod(U, x_A_x.Right()), adjoint(U));
   x_H_x.Right() = prod(prod(U, x_H_x.Right()), adjoint(U));
   x_Hb_x.Right() = prod(prod(U, x_Hb_x.Right()), adjoint(U));
   x_B_y.Right() = prod(U, x_B_y.Right());
   x_y.Right() = prod(U, x_y.Right(), Ident);

   this->DebugCheckBasis();
   return Info;
}

double FunctionalSolver::Solve(int MaxIterations)
{
   if (LinearSolver)
   {
      int Iter = MaxIterations;
      int m = MaxIterations+2;  // subspace size
      double Resid = 1e-10;

      MatrixOperator yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

      GmRes(x.Center(),
            SuperblockMultiply(conj(H.Center()),
                               x_H_x.Left(),
                               x_H_x.Right()),
            yprime,
            m, Iter, Resid,
            LinearAlgebra::Identity<MatrixOperator>());

      IterationNumMultiplies = Iter;
      // calculate the functional value
      std::complex<double> h = inner_prod(x.Center(), operator_prod(A.Center(), x_A_x.Left(), x.Center(), herm(x_A_x.Right())));
      std::complex<double> g = inner_prod(x.Center(), operator_prod(B.Center(), x_B_y.Left(), y.Center(), herm(x_B_y.Right())));
      DEBUG_TRACE(h)(g);
      double f = h.real() + 2.0 * g.real();
      return f;
   }
   else
   {
      if (SquareMeHarder)
      {
         MatrixOperator yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));

         double f = FunctionalMinimizeWithResidVector(x.Center(),
                                                      SuperblockMultiplySquare(conj(H.Center()),
                                                                               x_H_x.Left(),
                                                                               x_H_x.Right(),
                                                                               conj(Hb.Center()),
                                                                               x_Hb_x.Left(),
                                                                               x_Hb_x.Right()),
                                                      operator_prod(-1.0*conj(Hb.Center()),
                                                                    x_Hb_x.Left(),
                                                                    yprime,
                                                                    herm(x_Hb_x.Right())),
                                                      MaxIterations,
                                                      SuperblockMultiply(conj(H.Center()),
                                                                         x_H_x.Left(),
                                                                         x_H_x.Right()),
                                                      yprime,
                                                      norm_frob_sq(yprime),
                                                      Precision,
                                                      TwoStepKrylov,
                                                      MinIterations
                                                      );

         //         TRACE(f+norm_frob_sq(yprime));
         std::complex<double> h = inner_prod(x.Center(), operator_prod(A.Center(), x_A_x.Left(), x.Center(), herm(x_A_x.Right())));
         std::complex<double> g = inner_prod(x.Center(), operator_prod(B.Center(), x_B_y.Left(), y.Center(), herm(x_B_y.Right())));
         DEBUG_TRACE(h)(g);
         f = h.real() + 2.0 * g.real();

         IterationNumMultiplies = MaxIterations;
         return f;
      }
      else
      {
         if (UseResid)
         {
            MatrixOperator yprime = triple_prod(x_y.Left(), y.Center(), herm(x_y.Right()));
            double f = FunctionalMinimizeWithResidVector(x.Center(),
                                                         SuperblockMultiply(conj(A.Center()),
                                                                            x_A_x.Left(),
                                                                            x_A_x.Right()),
                                                         operator_prod(conj(B.Center()),
                                                                       x_B_y.Left(),
                                                                       y.Center(),
                                                                       herm(x_B_y.Right())),
                                                         MaxIterations,
                                                         SuperblockMultiply(conj(H.Center()),
                                                                            x_H_x.Left(),
                                                                            x_H_x.Right()),
                                                         yprime,
                                                         norm_frob_sq(yprime),
                                                         Precision,
                                                         TwoStepKrylov,
                                                         MinIterations
                                                         );

            IterationNumMultiplies = MaxIterations;
            return f;
         }
         else
         {
            double f = FunctionalMinimizeWithPre(x.Center(),
                                                 SuperblockMultiply(conj(A.Center()),
                                                                    x_A_x.Left(),
                                                                    x_A_x.Right()),
                                                 operator_prod(conj(B.Center()),
                                                               x_B_y.Left(),
                                                               y.Center(),
                                                               herm(x_B_y.Right())),
                                                 MaxIterations,
                                                 SuperblockMultiply(conj(H.Center()),
                                                                    x_H_x.Left(),
                                                                    x_H_x.Right()),
                                                 Precision,
                                                 TwoStepKrylov,
                                                 MinIterations
                                                 );

            IterationNumMultiplies = MaxIterations;
            return f;
         }
      }
   }
}

void FunctionalSolver::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(x_A_x.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_A_x.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_A_x.Left().Basis2(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_A_x.Right().Basis2(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_H_x.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_H_x.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_H_x.Left().Basis2(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_H_x.Right().Basis2(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_B_y.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_B_y.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_B_y.Left().Basis2(), y.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_B_y.Right().Basis2(), y.Center().Basis2());
}
