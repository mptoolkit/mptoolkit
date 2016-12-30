// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testregularize.cpp
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

#include "quantumnumbers/su2.h"
#include "quantumnumbers/u1.h"
#include "tensor/regularize.h"

using namespace QuantumNumbers;
using namespace Tensor;

int main()
{
   SymmetryList Q("S:SU(2)");
   QNConstructor<SU2> QN(Q);

   VectorBasis b1(Q);
   b1.push_back(QN(0.5), 5);
   b1.push_back(QN(0.5), 3);
   b1.push_back(QN(1), 4);
   b1.push_back(QN(0), 3);
   b1.push_back(QN(1), 1);
   b1.push_back(QN(0.5), 2);

   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> M = Regularize(b1);

   TRACE(M);
}
