
#include "blas/sparsematrix.h"

int main()
{
   blas::SparseMatrix<double> A(3,3);
   A.insert(2,2,3.4);

   std::cout << A << '\n';
}
