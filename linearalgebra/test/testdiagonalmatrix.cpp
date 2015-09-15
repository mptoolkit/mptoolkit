
#include "linearalgebra/diagonalmatrix.h"
#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   DiagonalMatrix<double> M(5,5);
   M.diagonal() = range(1,6);

   TRACE(M);

   Matrix<double> Mat = M;
   Mat(0,1) = 1;
   Mat(1,0) = 2;
   Mat(0,2) = 3;
   Mat(0,4) = 7;
   TRACE(Mat);

   Matrix<double> Mat2 = Mat*M;
   
   TRACE(Mat);
   //   TRACE(Matrix<double>(M*Mat));
   TRACE(Matrix<double>(Mat*M));

   TRACE(M*M);

   TRACE(norm_frob_sq(M*M));

   //   CHECK_EQUAL(norm_frob_sq(M*M), norm_frob_sq(Mat2));
}
