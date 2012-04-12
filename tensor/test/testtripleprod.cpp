// -*- C++ -*- $Id$

#include "tensor/tensor.h"
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

   IrredTensor<double> I(B, QN(0));
   I(0,0) = 1;
   I(1,1) = 1;
   I(2,2) = 1;

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

   CHECK_CLOSE(triple_prod(herm(N), I, P), scalar_prod(herm(N), P));

   CHECK_CLOSE(triple_prod(N, I, herm(P)), scalar_prod(N, herm(P)));
}
