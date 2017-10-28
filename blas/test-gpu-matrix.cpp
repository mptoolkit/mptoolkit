
#include "cuda/gpu_matrix.h"
#include "blas/matrix.h"
#include "cuda/cuda-setup.h"
#include "blas/vector.h"
#include "cuda/gpu_vector.h"

int main()
{
   cuda::setup_cuda();
   cublas::setup_cublas_thread();

   blas::Matrix<double> A({{1,2,3},
                           {4,5,6},
                           {7,8,9}});
   blas::Matrix<double> B({{10,11,12},
                           {13,14,15},
                           {16,17,18}});
   blas::Matrix<double> C(3,3);

   cublas::gpu_matrix<double> gA(3,3);
   cublas::gpu_matrix<double> gB(3,3);
   cublas::gpu_matrix<double> gC(3,3);

   // non-blocking set - cannot change A or B until operation completes
   set(gA, A);
   set(gB, B);

   gA.row(1) = 2*gB.row(0);

   std::cout << get_wait(gA) << '\n';
   std::cout << get_wait(gB) << '\n';

   gC = 2 * gA * gB;
   gC += 3 * gA * herm(gB);

   C = get_wait(gC);

   std::cout << C << '\n';

   blas::Vector<double> x({1.0, 2.0, 3.0});
   blas::Vector<double> y(3);

   cublas::gpu_vector<double> gx(3);
   cublas::gpu_vector<double> gy(3);

   set(gx, x);

   gy = 0.1*gC*gx;

   y = get_wait(gy);
   std::cout << y << '\n';

   y = get_wait(gA.diagonal());
   std::cout << y << '\n';

   cuda::gpu_ref<double> r = cuda::allocate_gpu_ref<double>();
   //   vector_sum(gA.diagonal(), r);
   trace(gA, r);

   std::cout << get_wait(r) << '\n';

   gC.row(0) = 2*gA*gy;

   std::cout << get_wait(gC) << '\n';

}
