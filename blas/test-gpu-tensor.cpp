
#include "cuda/gpu_matrix.h"
#include "blas/matrix.h"
#include "cuda/cuda-setup.h"
#include "blas/vector.h"
#include "cuda/gpu_vector.h"
#include "tensor/tensor.h"
#include "blas/matrix-blas.h"
#include <chrono>

using QuantumNumbers::QuantumNumber;
using Tensor::VectorBasis;

int main()
{
   std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

   cuda::setup_cuda();

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

   Tensor::IrredTensor<blas::gpu_matrix<double>, VectorBasis, VectorBasis>
      gM(B, B, QuantumNumber(SL));

   set(gM, M);

   Tensor::IrredTensor<blas::gpu_matrix<double>, VectorBasis, VectorBasis>
      gN(B, B, QuantumNumber(SL));

   set(gN, M);

#if 1
   Tensor::IrredTensor<blas::gpu_matrix<double>, VectorBasis, VectorBasis>
      gK = gN * gM;

   for (int i = 0; i < 1000; ++i)
   {
      if (i % 10 == 0)
         std::cout << i << '\n';
      add_prod(1.0, gM, gM, gK);
      //gK += gM*gM;
   }
#else
   Tensor::IrredTensor<blas::Matrix<double>, VectorBasis, VectorBasis>
      K = M * M;

   for (int i = 0; i < 1000; ++i)
   {
      if (i % 10 == 0)
         std::cout << i << '\n';
      add_prod(1.0, M,M, K);
      //K += M*M;
   }
#endif

   std::cout << "now\n";

   std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
   std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;

   cuda::device_synchronize();
   std::cout << "now\n";

   end = std::chrono::steady_clock::now();
   std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() <<std::endl;

}
