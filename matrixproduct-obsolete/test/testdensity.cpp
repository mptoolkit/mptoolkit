// -*- C++ -*- $Id$

#include "matrixproduct/density.h"
#include "quantumnumbers/su2.h"

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   BasisList B(Symmetry);
   B.push_back(QN(1));
   B.push_back(QN(1));

   SimpleOperator M(B, B, QN(0));
   M(0,0) = 0;
   M(0,1) = 0.5;
   M(1,0) = 0.5;
   M(1,1) = 0;

   DensityMatrix<SimpleOperator> DM(M);

   DM.DensityMatrixReport(std::cout);
}

