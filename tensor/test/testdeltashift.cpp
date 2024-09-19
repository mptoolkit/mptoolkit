// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/test/testdeltashift.cpp
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

#include "tensor/tensor.h"
#include "tensor/deltabasis.h"
#include "quantumnumbers/su2.h"
#include <iostream>
#include <iomanip>

using namespace Tensor;

using LinearAlgebra::equal;

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   BasisList B(Symmetry);
   B.push_back(QN(0));
   B.push_back(QN(1));
   B.push_back(QN(2));

   IrredTensor<double> N(B, QN(1));
   N(0,1) = 2;
   N(1,2) = 4;
   N(2,1) = 6;

   IrredTensor<double> P(B, QN(1));
   P(0,1) = 7;
   P(1,0) = 3;
   P(2,1) = 5;

   TRACE(inner_prod(N, N));
   TRACE(inner_prod(P, P));

   TRACE(N);
   TRACE(N.Basis1())(N.Basis2());

   double nNorm = norm_frob(N);

   IrredTensor<double> Ns = delta_shift(N, QuantumNumbers::Projection(N.GetSymmetryList(), "1"));

   TRACE(Ns);
   TRACE(Ns.Basis1())(Ns.Basis2());
   TRACE(inner_prod(Ns, Ns));

   CHECK_EQUAL(norm_frob(Ns), nNorm);

}
