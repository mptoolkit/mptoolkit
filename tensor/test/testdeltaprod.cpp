// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/test/testdeltaprod.cpp
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
// $Id$

// this is a test of a case where the tensor product is giving unexpected results.  Probably
// not a problem with tensor_prod, but more likely some PBKAC.

#include "tensor/tensorproduct.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"

using namespace Tensor;
using namespace QuantumNumbers;

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   // helper function object for creating quantum numbers at compile time
   QNConstructor<SU2> QN(Symmetry);

   BasisList TestBasis1(Symmetry);
   TestBasis1.push_back(QN(15));

   BasisList TestBasis2(Symmetry);
   TestBasis2.push_back(QN(14.5));
   TestBasis2.push_back(QN(15.5));

   IrredTensor<double> A(TestBasis1, TestBasis2, QN(0.5));

   //A(0,0) = 0.866025;
   A(0,0) = 1;
   A(0,1) = 100;
   //A(0,1) = 1.11803;
   //A(2,9) = 1.08012;
   //A(3,6) = 0.912871;
   //A(4,7) = 1.09545;
   //A(5,8) = 1.09545;

   TRACE(scalar_prod(herm(A), A));
   TRACE(scalar_prod(A, herm(A)));

   BasisList DeltaBasis(Symmetry);
   DeltaBasis.push_back(QN(15));
   IrredTensor<double> Delta(DeltaBasis, DeltaBasis, QN(0));
   Delta(0,0) = 1;

   IrredTensor<double> Res = Tensor::tensor_prod(A, Delta, A.TransformsAs());
   TRACE(scalar_prod(herm(Res), Res))(scalar_prod(Res, herm(Res)));

   TRACE(Res.Basis1())(Res.Basis2())(Res);
}
