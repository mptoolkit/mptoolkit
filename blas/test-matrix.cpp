
#include "blas/matrix.h"
#include "blas/matrix-blas.h"

int main()
{
   blas::Matrix<double> A({{1,2,3},{4,5,6},{7,8,9}});
   blas::Matrix<double> B({{10,11,12},{13,14,15},{16,17,18}});
   blas::Matrix<double> C(3,3);

   C = 2 * A * B;
   C += 3 * A * herm(B);

   std::cout << C << '\n';
}
