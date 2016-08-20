// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/test/testtensor2.cpp
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

   TRACE(B);

   IrredTensor<double> M(B, QN(0));
   M(0,0) = 1;
   M(1,1) = 2;
   M(2,2) = 3;

   TRACE(M)(norm_2(M));

   IrredTensor<double> N(B, QN(1));
   N(0,1) = 2;
   N(1,2) = 4;
   N(2,1) = 6;

   TRACE(N)(norm_2(N));

   IrredTensor<double> P(B, QN(1));
   P(1,0) = 3;
   P(2,1) = 5;
   P(0,1) = 7;

   TRACE(P)(norm_2(P));

   std::cout << scalar_prod(N, P) << std::endl;

   TRACE(prod(adjoint(N), P, QN(0)));
   TRACE(prod(N, adjoint(P), QN(0)));
   TRACE(scalar_product(herm(N), P));
   TRACE(scalar_product(N, herm(P)));

   std::cout << prod(M, N) << std::endl;
   std::cout << prod(M, N, QN(1), TestFunc()) << std::endl;

   show_projections(std::cout, prod(M, N, QN(1), TestFunc()));

   std::cout << tensor_prod(M, N, QN(1)) << std::endl;

   std::cout << tensor_sum(N, P) << std::endl;
   std::cout << tensor_row_sum(N, P) << std::endl;
   std::cout << tensor_col_sum(N, P) << std::endl;

}
