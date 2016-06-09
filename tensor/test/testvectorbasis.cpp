// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testvectorbasis.cpp
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

#include "tensor/basis.h"
#include "tensor/tensorsum.h"
#include "tensor/tensorproduct.h"
#include "quantumnumbers/u1.h"
#include "quantumnumbers/su2.h"

using namespace Tensor;
using namespace QuantumNumbers;

int main()
{
   SymmetryList Symmetry("N:U(1),S:SU(2)");
   QNConstructor<U1,SU2> QN(Symmetry);

   BasisList B1(Symmetry);
   B1.push_back(QN(0, 0));   // empty
   B1.push_back(QN(1, 0.5)); // single
   B1.push_back(QN(2, 0));   // double
   B1.push_back(QN(0, 0));   // another empty

   int Dims[4] = {2,2,1,1};

   VectorBasis VB1(B1, &Dims[0], &Dims[4]);

   CHECK_EQUAL(VB1, VB1);
   CHECK_EQUAL(VB1.dim(0), 2);
   CHECK_EQUAL(VB1.dim(1), 2);
   CHECK_EQUAL(VB1.dim(2), 1);
   CHECK_EQUAL(VB1.dim(3), 1);
   CHECK_EQUAL(VB1.total_dimension(), 6);
   CHECK_EQUAL(VB1.total_degree(), 8);

   BasisList B2(Symmetry);
   B2.push_back(QN(1, 0.5));   
   B2.push_back(QN(2, 0.5));
   int Dims2[2] = {2,1};

   VectorBasis VB2(B2, &Dims2[0], &Dims2[2]);

   CHECK_EQUAL(VB2, VB2);

   CHECK(VB1 != VB2);

   SumBasis<VectorBasis> VSum(VB1, VB2);
   CHECK_EQUAL(VSum.Basis().total_dimension(), VB1.total_dimension()+VB2.total_dimension());
   CHECK_EQUAL(VSum.Basis().total_degree(), VB1.total_degree()+VB2.total_degree());

   ProductBasis<VectorBasis, VectorBasis> VProd(VB1, VB2);
   CHECK_EQUAL(VProd.Basis().total_degree(), VB1.total_degree()*VB2.total_degree());
}
