// -*- C++ -*- $Id$

#include "matrixproduct/periodicwavefunction.h"
#include "quantumnumbers/su2.h"

int main()
{
   // construct an AKLT state
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);

   BasisList LocalBasis(Symmetry);
   LocalBasis.push_back(QN(1));  // S=1 chain

   VectorBasis MatrixBasis(Symmetry);
   MatrixBasis.push_back(QN(0.5));    // spin 1/2 edge states

   PeriodicWavefunction AKLT(100, LocalBasis, MatrixBasis);

   AKLT.Q() = MatrixOperator::make_identity(MatrixBasis);
   AKLT[0] = MatrixOperator(MatrixBasis, QN(1));
   AKLT[0](0,0) = LinearAlgebra::Matrix<double>(1,1);
   AKLT[0](0,0)(0,0) = sqrt(2.0); // sqrt(s*(s+1)) = sqrt(2.0)

   TRACE(AKLT);
}
