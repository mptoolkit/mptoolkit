// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ddmrg.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "ddmrg.h"
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

Solver::Solver(CenterWavefunction const& Psi_, SplitOperator const& A_, CenterWavefunction const& Rhs_, double Freq_, double Broad_)
   : x(Psi_), A(A_), y(Rhs_), Frequency(Freq_), Broadening(Broad_), Ident(Psi_.GetSymmetryList())
{
   // construct the OpMatrix elements for the right hand side

   typedef MPStateComponent Component;

   // Shift A so the Center matrix lines up with that of x
   while (A.LeftSize() > x.LeftSize()) A.RotateLeft();
   while (A.RightSize() > x.RightSize()) A.RotateRight();

   // Shift y so the Center matrix lines up with that of x
   while (y.LeftSize() > y.LeftSize()) y.RotateLeft();
   while (y.RightSize() > y.RightSize()) y.RotateRight();


   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());

   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);

   x_A_x.PushRight(make_vacuum_state(x.GetSymmetryList()));
   x_y.PushRight(VacIdent);

   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      x_A_x.PushRight(operator_prod(A.LookupRight(Loc),
                                    x.LookupRight(Loc),
                                    x_A_x.Right(),
                                    herm(x.LookupRight(Loc))));

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
   x_y.PushLeft(VacIdent);

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      x_A_x.PushLeft(operator_prod(herm(A.LookupLeft(Loc)),
                                   herm(x.LookupLeft(Loc)),
                                   x_A_x.Left(), 
                                   x.LookupLeft(Loc)));

      x_y.PushLeft(operator_prod(herm(x.LookupLeft(Loc)),
                                 x_y.Left(), 
                                 y.LookupLeft(Loc)));
   }

   this->DebugCheckBasis();
}

Solver::~Solver()
{
}

void Solver::ExpandLeft()
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

   x_y.PopLeft();
   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));

   this->DebugCheckBasis();
}

void Solver::ExpandRight()
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

   x_y.PopRight();
   x_y.PushRight(operator_prod(x.Right(), x_y.Right(), herm(y.Right())));

   this->DebugCheckBasis();
}

void Solver::ShiftRightAndExpand()
{
   // Rotate the wavefunctions to the right.  A and y are fixed and do not
   // need the basis expanded

   A.RotateRight();
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

   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));
   x_y.PopRight();

   this->DebugCheckBasis();
}

void Solver::ShiftLeftAndExpand()
{
   // Rotate the wavefunctions to the left.  A and y are fixed and do not
   // need the basis expanded
   A.RotateLeft();
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

   x_y.PushRight(operator_prod(x.Right(),
                               x_y.Right(),
                               herm(y.Right())));
   x_y.PopLeft();

   this->DebugCheckBasis();
}

TruncationInfo Solver::TruncateLeft(int MinStates, int MaxStates, double MinTrunc, double CFactor)
{
   // Calculate the density matrix for x
   MatrixOperator Rho_x = scalar_prod(x.Center(), herm(x.Center()));

   MatrixOperator Rho = Rho_x * (1.0 / trace(Rho_x));

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
      = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                           MinStates,
                                           MaxStates,
                                           MinTrunc,
                                           Info);

   // truncate
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   x.Center() = prod(U, x.Center(), Ident);
   x.Left() = prod(x.Left(), adjoint(U));

   x_A_x.Left() = prod(prod(U, x_A_x.Left()), adjoint(U));
   x_y.Left() = prod(U, x_y.Left(), Ident);

   this->DebugCheckBasis();
   return Info;
}

TruncationInfo Solver::TruncateRight(int MinStates, int MaxStates, double MinTrunc, double CFactor)
{
   MatrixOperator Rho_x = scalar_prod(herm(x.Center()), x.Center());

   MatrixOperator Rho = Rho_x * (1.0 / trace(Rho_x));

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
      = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                           MinStates,
                                           MaxStates,
                                           MinTrunc,
                                           Info);
   // truncate
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), E);

   x.Center() = prod(x.Center(), adjoint(U), Ident);
   x.Right() = prod(U, x.Right());

   x_A_x.Right() = prod(prod(U, x_A_x.Right()), adjoint(U));
   x_y.Right() = prod(U, x_y.Right(), Ident);

   this->DebugCheckBasis();
   return Info;
}

void Solver::DebugCheckBasis() const
{
   DEBUG_CHECK_EQUAL(x_A_x.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_A_x.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_A_x.Left().Basis2(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_A_x.Right().Basis2(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_y.Left().Basis1(), x.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_y.Right().Basis1(), x.Center().Basis2());

   DEBUG_CHECK_EQUAL(x_y.Left().Basis2(), y.Center().Basis1());
   DEBUG_CHECK_EQUAL(x_y.Right().Basis2(), y.Center().Basis2());
}
