// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/wavefunc-utils.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "wavefunc-utils.h"
#include "linearalgebra/matrix_utility.h"

//
// TransformStack
//

void InitializeTransformStack(TransformOperator& x_y,
                              CenterWavefunction const& x,
                              CenterWavefunction const& y)
{
   x_y = TransformOperator();
   // make sure the two input wavefunctions are at the same rotation
   CHECK_EQUAL(x.LeftSize(), y.LeftSize());
   CHECK_EQUAL(x.RightSize(), y.RightSize());
   // make sure the two states transform the same way
   CHECK_EQUAL(x.LookupLeft(0).Basis1(), y.LookupLeft(0).Basis1());
   // Make sure the state is irreducible.  This probably isn't strictly necessary?
   CHECK_EQUAL(x.LookupLeft(0).Basis1().size(), 1);

   QuantumNumber Ident(x.GetSymmetryList());

   // The initial matrix elements represent the vacuum basis on the right hand side
   BasisList Vacuum = make_vacuum_basis(x.GetSymmetryList());

   MatrixOperator VacIdent = MatrixOperator(VectorBasis(Vacuum), VectorBasis(Vacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   x_y.PushRight(VacIdent);

   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      x_y.PushRight(operator_prod(x.LookupRight(Loc),
                                  x_y.Right(),
                                  herm(y.LookupRight(Loc))));
   }

   //Now start from the left hand side and work inwards
   BasisList LeftVacuum(x.GetSymmetryList());
   LeftVacuum.push_back(x.LookupLeft(0).Basis1()[0]);

   VacIdent = MatrixOperator(VectorBasis(LeftVacuum), VectorBasis(LeftVacuum), Ident);
   VacIdent(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   x_y.PushLeft(VacIdent);

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      x_y.PushLeft(operator_prod(herm(x.LookupLeft(Loc)),
                                 x_y.Left(),
                                 y.LookupLeft(Loc)));
   }
}

void TransformStackRotateLeft(TransformOperator& x_y,
                              CenterWavefunction const& x,
                              CenterWavefunction const& y)
{
   x_y.PushRight(operator_prod(x.Right(), x_y.Right(), herm(y.Right())));
   x_y.PopLeft();
}

void TransformStackRotateRight(TransformOperator& x_y,
                               CenterWavefunction const& x,
                               CenterWavefunction const& y)
{
   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));
   x_y.PopRight();
}

void TransformStackUpdateLeft(TransformOperator& x_y,
                              CenterWavefunction const& x,
                              CenterWavefunction const& y)
{
   x_y.PopLeft();
   x_y.PushLeft(operator_prod(herm(x.Left()), x_y.Left(), y.Left()));
}

void TransformStackUpdateRight(TransformOperator& x_y,
                               CenterWavefunction const& x,
                               CenterWavefunction const& y)
{
   x_y.PopRight();
   x_y.PushRight(operator_prod(x.Right(), x_y.Right(), herm(y.Right())));
}

//
// SuperblockStack
//

void InitializeSuperblockStack(SuperblockOperator& x_A_y,
                               CenterWavefunction const& x,
                               SplitOperator& A,
                               CenterWavefunction const& y)
{
   x_A_y = SuperblockOperator();

   while (A.LeftSize() < x.LeftSize())
      A.RotateRight();
   while (A.LeftSize() > x.LeftSize())
      A.RotateLeft();

   // The initial matrix elements represent the vacuum basis on the right hand side
   x_A_y.PushRight(make_vacuum_state(x.GetSymmetryList()));

   // apply the Right matrices
   for (int Loc = 0; Loc < x.RightSize(); ++Loc)
   {
      x_A_y.PushRight(operator_prod(A.LookupRight(Loc),
                                    x.LookupRight(Loc),
                                    x_A_y.Right(),
                                    herm(y.LookupRight(Loc))));
   }

   //Now start from the left hand side and work inwards
   x_A_y.PushLeft(make_vacuum_state(x.LookupLeft(0).Basis1()[0]));

   for (int Loc = 0; Loc < x.LeftSize(); ++Loc)
   {
      x_A_y.PushLeft(operator_prod(herm(A.LookupLeft(Loc)),
                                   herm(x.LookupLeft(Loc)),
                                   x_A_y.Left(),
                                   y.LookupLeft(Loc)));
   }
}

void SuperblockStackRotateLeft(SuperblockOperator& x_A_y,
                               CenterWavefunction const& x,
                               SplitOperator const& A,
                               CenterWavefunction const& y)
{
   x_A_y.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_y.Right(),
                                 herm(y.Right())));
   x_A_y.PopLeft();
}

void SuperblockStackRotateRight(SuperblockOperator& x_A_y,
                                CenterWavefunction const& x,
                                SplitOperator const& A,
                                CenterWavefunction const& y)
{
   x_A_y.PushLeft(operator_prod(herm(A.Left()),
                                herm(x.Left()),
                                x_A_y.Left(),
                                y.Left()));
   x_A_y.PopRight();
}

void SuperblockStackUpdateLeft(SuperblockOperator& x_A_y,
                               CenterWavefunction const& x,
                               SplitOperator const& A,
                               CenterWavefunction const& y)
{
   x_A_y.PopLeft();
   x_A_y.PushLeft(operator_prod(herm(A.Left()),
                                herm(x.Left()),
                                x_A_y.Left(),
                                y.Left()));
}

void SuperblockStackUpdateRight(SuperblockOperator& x_A_y,
                                CenterWavefunction const& x,
                                SplitOperator const& A,
                                CenterWavefunction const& y)
{
   x_A_y.PopRight();
   x_A_y.PushRight(operator_prod(A.Right(),
                                 x.Right(),
                                 x_A_y.Right(),
                                 herm(y.Right())));
}
