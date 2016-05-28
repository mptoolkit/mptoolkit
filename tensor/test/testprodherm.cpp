// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testprodherm.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "tensor/tensor.h"
#include "tensor/tensorproduct.h"
#include "tensor/tensorsum.h"
#include "quantumnumbers/su2.h"
#include <iostream>
#include <iomanip>

using namespace Tensor;

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

   IrredTensor<double> N(B, QN(1));
   N(0,1) = 2;
   N(1,2) = 4;
   N(2,1) = 6;

   IrredTensor<double> P(B, QN(1));
   P(1,0) = 3;
   P(2,1) = 5;
   P(0,1) = 7;

   CHECK_CLOSE(N*herm(M), N*adjoint(M));
   CHECK_CLOSE(P*herm(M), P*adjoint(M));

   CHECK_CLOSE(herm(M)*N, adjoint(M)*N);
   CHECK_CLOSE(herm(M)*P, adjoint(M)*P);
}
