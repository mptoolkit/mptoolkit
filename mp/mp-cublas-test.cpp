
#include "cuda/cuda-setup.h"
#include "cuda/cublas.h"

int main()
{
   cuda::setup_cuda();

   TRACE("start");
   cublas::handle H = cublas::handle::create();


   TRACE("arena");
   cuda::arena Arena = cuda::get_block_allocator();

   LinearAlgebra::Matrix<double> M(3,3,0.0);
   M(0,0) = 1;
   M(0,1) = 2;
   M(0,2) = 3;
   M(1,0) = 4;
   M(1,1) = 5;
   M(1,2) = 6;
   M(2,0) = 7;
   M(2,1) = 8;
   M(2,2) = 9;

   TRACE("matrix1");
   cublas::gpu_matrix<double> A(3,3,Arena);
   TRACE("matrix2");
   cublas::gpu_matrix<double> B(3,3,Arena);

   TRACE("set1");
   cublas::set_wait(A, M);

   LinearAlgebra::Matrix<double> AA = get_wait(A);
   std::cout << AA << '\n';


   M *= 2;

   TRACE("set1");
   cublas::set_wait(B, M);

   TRACE("matrix3");
   cublas::gpu_matrix<double> C(3,3,Arena);

   TRACE("gemm");
   cublas::gemm(H, 1.0, A, B, 0.0, C);

   LinearAlgebra::Matrix<double> N = get_wait(C);

   std::cout << N << '\n';
}
