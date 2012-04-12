
#include "matrixproduct/mpstate.h"
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

   MPStateComponent AKLT(LocalBasis, MatrixBasis, MatrixBasis);
   AKLT[0] = MatrixOperator(MatrixBasis, QN(1));
   AKLT[0](0,0) = LinearAlgebra::Matrix<double>(1,1);
   AKLT[0](0,0)(0,0) = sqrt(2.0); // sqrt(s*(s+1)) = sqrt(2.0)

   TRACE(LinearAlgebra::scalar_direct_prod(herm(AKLT), AKLT));
   TRACE(LinearAlgebra::scalar_direct_prod(AKLT, herm(AKLT)));
}
