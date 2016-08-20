// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testtensorsum.cpp
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

   IrredTensor<double> N(B, QN(1));
   N(0,1) = 2;
   N(1,2) = 4;
   N(2,1) = 6;

   IrredTensor<double> P(B, QN(1));
   P(0,1) = 7;
   P(1,0) = 3;
   P(2,1) = 5;

   SumBasis<BasisList> SB(B,B);
   IrredTensor<double> NP(SB, QN(1));
   NP(0,1) = 2;
   NP(1,2) = 4;
   NP(2,1) = 6;
   NP(3,4) = 7;
   NP(4,3) = 3;
   NP(5,4) = 5;

   std::vector<IrredTensor<double> > NPacc(2);
   NPacc[0] = N;
   NPacc[1] = P;


   CHECK_EQUAL(SB, tensor_sum(N, P).Basis1());
   CHECK_EQUAL(SB, tensor_sum(N, P).Basis2());
   CHECK_CLOSE(NP, tensor_sum(N, P));
   CHECK_CLOSE(NP, tensor_sum(N, P, SB));

   CHECK_EQUAL(SB, tensor_accumulate(NPacc.begin(), NPacc.end(), SB, SB).Basis1());
   CHECK_EQUAL(SB, tensor_accumulate(NPacc.begin(), NPacc.end(), SB, SB).Basis2());
   CHECK_CLOSE(NP, tensor_accumulate(NPacc.begin(), NPacc.end(), SB, SB));

   IrredTensor<double> NProw(B, SB, QN(1));
   NProw(0,1) = 2;
   NProw(1,2) = 4;
   NProw(2,1) = 6;
   NProw(0,4) = 7;
   NProw(1,3) = 3;
   NProw(2,4) = 5;

   CHECK_EQUAL(B, tensor_row_sum(N, P).Basis1());
   CHECK_EQUAL(SB, tensor_row_sum(N, P).Basis2());
   CHECK_CLOSE(NProw, tensor_row_sum(N, P));
   CHECK_CLOSE(NProw, tensor_row_sum(N, P, SB));

   CHECK_EQUAL(B, tensor_row_accumulate(NPacc.begin(), NPacc.end(), SB).Basis1());
   CHECK_EQUAL(SB, tensor_row_accumulate(NPacc.begin(), NPacc.end(), SB).Basis2());
   CHECK_CLOSE(NProw, tensor_row_accumulate(NPacc.begin(), NPacc.end(), SB));

   IrredTensor<double> NPcol(SB, B, QN(1));
   NPcol(0,1) = 2;
   NPcol(1,2) = 4;
   NPcol(2,1) = 6;
   NPcol(3,1) = 7;
   NPcol(4,0) = 3;
   NPcol(5,1) = 5;

   CHECK_EQUAL(SB, tensor_col_sum(N, P).Basis1());
   CHECK_EQUAL(B, tensor_col_sum(N, P).Basis2());
   CHECK_CLOSE(NPcol, tensor_col_sum(N, P));
   CHECK_CLOSE(NPcol, tensor_col_sum(N, P, SB));

   CHECK_EQUAL(SB, tensor_col_accumulate(NPacc.begin(), NPacc.end(), SB).Basis1());
   CHECK_EQUAL(B, tensor_col_accumulate(NPacc.begin(), NPacc.end(), SB).Basis2());
   CHECK_CLOSE(NPcol, tensor_col_accumulate(NPacc.begin(), NPacc.end(), SB));
}
