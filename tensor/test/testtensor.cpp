// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testtensor.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "tensor/tensor.h"
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
   CHECK_EQUAL(N.size1(), 3);
   CHECK_EQUAL(N.size2(), 3);
   CHECK_EQUAL(N.TransformsAs(), QN(1));
   CHECK_EQUAL(size1(N.data()), 3);
   CHECK_EQUAL(size2(N.data()), 3);
   N(0,1) = 2;
   N(1,2) = 4;
   N(2,1) = 6;

   IrredTensor<double> P(B, QN(1));
   P(0,1) = 7;
   P(1,0) = 3;
   P(2,1) = 5;

   IrredTensor<double> NpP = N+P;
   CHECK_EQUAL(NpP.TransformsAs(), QN(1));

   CHECK_EQUAL(NpP(0,1), 9);
   CHECK_EQUAL(NpP(1,0), 3);
   CHECK_EQUAL(NpP(1,2), 4);
   CHECK_EQUAL(NpP(2,1), 11);

   // make sure output works
   TRACE(NpP);
   TRACE(show_projections(NpP));

   IrredTensor<double> Nc = - (P * 2.0 -  2.0 * NpP);
   N *= 2.0;
   CHECK_CLOSE(Nc, N);

   Nc = transform(N, LinearAlgebra::Negate<double>());
   CHECK_CLOSE(Nc, -N);

   // trace
   IrredTensor<double> Tt(B, QN(0));
   Tt(0,0) = 2;
   Tt(1,1) = 10;
   Tt(2,2) = 100;
   CHECK_EQUAL(trace(Tt), 532);

   // norm_frob
   CHECK_CLOSE(norm_frob_sq(P), 201);
   CHECK_CLOSE(norm_frob_sq(Tt), 50304);
   CHECK_CLOSE(norm_frob_sq(N), 928);
   CHECK_CLOSE(norm_frob(N), 30.46309242345563);

   CHECK_CLOSE(inner_prod(P, P), norm_frob_sq(P));
   CHECK_CLOSE(inner_prod(P, N), 328);
   CHECK_CLOSE(inner_prod(N, P), 328);
}
