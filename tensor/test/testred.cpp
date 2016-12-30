// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testred.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "tensor/reducible.h"
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

   ReducibleTensor<double> N(B);
   CHECK_EQUAL(N.size1(), 3);
   CHECK_EQUAL(N.size2(), 3);
   CHECK_EQUAL(size1(N.project(QN(0)).data()), 3);
   CHECK_EQUAL(size2(N.project(QN(0)).data()), 3);
   N.project(QN(1))(0,1) = 2;
   N.project(QN(1))(1,2) = 4;
   N.project(QN(1))(2,1) = 6;

   TRACE(N);

   std::cout << show_projections(N) << '\n';
}
