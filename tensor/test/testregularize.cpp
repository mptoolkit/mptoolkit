// -*- C++ -*- $Id$

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
