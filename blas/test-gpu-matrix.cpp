
#include "cuda/gpu_matrix.h"
#include "blas/matrix.h"

int main()
{
   blas::Matrix<double> A({{1,2,3},{4,5,6},{7,8,9}});
   blas::Matrix<double> B({{10,11,12},{13,14,15},{16,17,18}});
   blas::Matrix<double> C(3,3);

   arena Mem = cuda::get_block_allocator();

   cublas::gpu_matrix<double> gA(3,3, Mem);
   cublas::gpu_matrix<double> gB(3,3, Mem);
   cublas::gpu_matrix<double> gC(3,3, Mem);

   gA = A;
   gB = B;

   gC = 2 * gA * gB;

   C = get_wait(gC);

   std::cout << C << '\n';
}
