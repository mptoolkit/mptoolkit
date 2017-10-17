
#include "cuda/gpu_matrix.h"
#include "blas/matrix.h"
#include "cuda/cuda-setup.h"
#include "blas/vector.h"
#include "cuda/gpu_vector.h"
#include "tensor/tensor.h"
#include "blas/matrix-blas.h"

using QuantumNumbers::QuantumNumber;
using Tensor::VectorBasis;

int main()
{
   cuda::setup_cuda();
   cublas::setup_cublas_thread();

   QuantumNumbers::SymmetryList SL("U:U(1)");
   Tensor::VectorBasis B(SL);
   B.push_back(QuantumNumber(SL, "0"), 1000);
   B.push_back(QuantumNumber(SL, "1"), 500);
   B.push_back(QuantumNumber(SL, "2"), 1000);

   Tensor::IrredTensor<blas::Matrix<double>, VectorBasis, VectorBasis>
      M(B, B, QuantumNumber(SL));

   M.insert(0, 0, blas::random_matrix<double>(1000,1000));
   M.insert(1, 1, blas::random_matrix<double>(500,500));
   M.insert(2, 2, blas::random_matrix<double>(1000,1000));

   std::cout << "now\n";

   Tensor::IrredTensor<cublas::gpu_matrix<double>, VectorBasis, VectorBasis>
      gM(B, B, QuantumNumber(SL));

   set(gM, M);

   Tensor::IrredTensor<cublas::gpu_matrix<double>, VectorBasis, VectorBasis>
      gN(B, B, QuantumNumber(SL));

   set(gN, M);

#if 1
   Tensor::IrredTensor<cublas::gpu_matrix<double>, VectorBasis, VectorBasis>
      gK = gN * gM;

   for (int i = 0; i < 1000; ++i)
   {
      if (i % 10 == 0)
         std::cout << i << '\n';
      add_prod(gM, gM, 1.0, gK);
      //gK += gM*gM;
   }
#else
   Tensor::IrredTensor<blas::Matrix<double>, VectorBasis, VectorBasis>
      K = M * M;

   for (int i = 0; i < 1000; ++i)
   {
      add_prod(M,M, 1.0, K);
      //K += M*M;
   }
#endif

   std::cout << "now\n";


   cuda::device_synchronize();
   std::cout << "now\n";

}
