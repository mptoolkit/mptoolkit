
#include "blas/matrix.h"
#include "blas/vector.h"
#include "blas/matrix-blas.h"

int main()
{
   blas::Matrix<double> A({{1,2,3},
                           {4,5,6},
                           {7,8,9}});
   blas::Matrix<double> B({{10,11,12},
                           {13,14,15},
                           {16,17,18}});
   blas::Matrix<double> C(3,3);

   A.row(0) = 2 * A.row(0);

   C = 2 * A * B;
   C += 3 * A * herm(B);

   std::cout << C << '\n';


   blas::Vector<double> x({1.0, 2.0, 3.0});
   blas::Vector<double> y(3);

   y = 2*C*x;

   std::cout << y << '\n';

   y = A.diagonal();
   std::cout << y << '\n';

}
