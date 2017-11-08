
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

   A.row(1) = 2 * B.row(0);

   std::cout << A << '\n';
   std::cout << B << '\n';

   C = 2 * A * B;
   C += 3 * A * herm(B);

   std::cout << C << '\n';

   blas::Vector<double> x({1.0, 2.0, 3.0});
   blas::Vector<double> y(3);

   y = 0.1*C*x;

   std::cout << y << '\n';

   std::cout << A.diagonal() << '\n';

   double r;
   trace(A, r);

   std::cout << r << '\n';

   C.row(0) = 2*A*y;

   std::cout << C << '\n';

   std::cout << "Testing lapack\n";
   C = herm(A);
   A += C;
   std::cout << A << '\n';
   blas::DiagonalizeSymmetric(A, x);

   std::cout << A << '\n' << x << '\n';

   C.clear();
   C.diagonal() = x;
   B = A * C;
   C = B * herm(A);

   // should be original matrix
   std::cout << C << '\n';

}
