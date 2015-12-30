
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/diagonalmatrix.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

int main()
{
   for (int i = 1; i <= 10; ++i)
   {
      for (int j = 1; j <= 10; ++j)
      {
	 Matrix<std::complex<double>> M = random_matrix<std::complex<double>>(i,j);

	 Matrix<std::complex<double>> U, Vh;
	 DiagonalMatrix<double> D;
	 
	 SingularValueDecompositionFull(M, U, D, Vh);
	 
	 CHECK(norm_frob(U*D*Vh - M) < 1E-10)(M)(U)(D)(Vh);
      }
   }
}
