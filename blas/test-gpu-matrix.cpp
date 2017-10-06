
#include "cuda/gpu_matrix.h"
#include "blas/matrix.h"
#include "cuda/cuda-setup.h"

int main()
{
   cuda::setup_cuda();
   cublas::setup_cublas_thread();

   blas::Matrix<double> A({{1,2,3},{4,5,6},{7,8,9}});
   blas::Matrix<double> B({{10,11,12},{13,14,15},{16,17,18}});
   blas::Matrix<double> C(3,3);

   cublas::gpu_matrix<double> gA(3,3);
   cublas::gpu_matrix<double> gB(3,3);
   cublas::gpu_matrix<double> gC(3,3);

   // non-blocking set - cannot change A or B until operation completes
   set(gA, A);
   set(gB, B);

   gC = 2 * gA * gB;
   gC += 3 * gA * herm(gB);

   C = get_wait(gC);

   std::cout << C << '\n';
}
