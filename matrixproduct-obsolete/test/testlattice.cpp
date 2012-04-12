// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "models/su2spinchains.h"

int main()
{
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);

   Lattice L;
   L.Append("1", CreateSU2SpinSite(0.5));
   L.Append("2", CreateSU2SpinSite(1));
   L.Append("3", CreateSU2SpinSite(1));
   L.Append("4", CreateSU2SpinSite(0.5));
}
