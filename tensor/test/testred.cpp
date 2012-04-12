// -*- C++ -*- $Id$

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
