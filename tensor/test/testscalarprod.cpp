// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testscalarprod.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "tensor/tensor.h"
#include "tensor/tensorproduct.h"
#include "tensor/tensorsum.h"
#include "quantumnumbers/su2.h"
#include <iostream>
#include <iomanip>

using namespace Tensor;

struct TestFunc
{
   typedef double first_argument_type;
   typedef double second_argument_type;
   typedef double result_type;
   typedef double value_type;

   result_type operator()(double x, double y) const { return x*y * 100; }
};

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   BasisList B(Symmetry);
   B.push_back(QN(0));
   B.push_back(QN(1));
   B.push_back(QN(2));

   IrredTensor<double> M(B, QN(0));
   M(0,0) = 1;
   M(1,1) = 2;
   M(2,2) = 3;

   CHECK_EQUAL(trace(M), 22);
   CHECK_EQUAL(norm_frob_sq(M), 58);

   IrredTensor<double> N(B, QN(1));
   N(0,1) = 2;
   N(1,2) = 4;
   N(2,1) = 6;

   IrredTensor<double> P(B, QN(1));
   P(1,0) = 3;
   P(2,1) = 5;
   P(0,1) = 7;

   IrredTensor<double> aNP = prod(adjoint(N), P, QN(0));
   IrredTensor<double> hNP = scalar_prod(herm(N), P);
   CHECK(equal(sqrt(3.0) * aNP, hNP))(sqrt(3.0) * aNP)(hNP);
   CHECK(equal(hNP, adjoint(scalar_prod(herm(P), N))));

   IrredTensor<double> NaP = prod(N, adjoint(P), QN(0));
   IrredTensor<double> NhP = scalar_prod(N, herm(P));
   CHECK(equal(sqrt(3.0) * NaP, NhP))(sqrt(3.0) * NaP)(NhP);
   CHECK(equal(NhP, adjoint(scalar_prod(P, herm(N)))));

   // test custom functor argument for the prod() function
   CHECK(equal(100 * prod(M, N, QN(1)), prod(M, N, QN(1), TestFunc()), 1.0e-10))
      (100 * prod(M, N, QN(1)))
      (prod(M, N, QN(1), TestFunc()));

}
