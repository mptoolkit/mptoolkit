// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cublas.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "cublas.h"

namespace blas
{

arena gpu_default_arena = blas::arena(new cuda::BlockAllocator(cuda::DefaultBlockMultiple, false));

} // namespace blas

namespace cublas
{

char const*
GetErrorName(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}

char const*
GetErrorString(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "cuBLAS: SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "cuBLAS: library not initialized";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "cuBLAS: memory allocation failed";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "cuBLAS: invalid value or parameter";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "cuBLAS: feature not supported by this architecture (possible double-precision?)";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "cuBLAS: invalid GPU memory mapping";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "cuBLAS: GPU kernel execution failed";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "cuBLAS: internal error";
    }

    return "<cublas-unknown>";
}

} // namespace cubas

namespace cuda
{

void
matrix_inner_prod(char Atrans, char Btrans, int M, int N,
		      cuda::const_gpu_ptr<double> A, int ldA,
		      cuda::const_gpu_ptr<double> B, int ldB,
		      cuda::gpu_ref<double>& r)
{
   if (Atrans == 'N' && Btrans == 'N')
   {
      cuda::gpu_buffer<double> Acc = cuda::allocate_gpu_temporary<double>(N);

      for (int i = 0; i < N; ++i)
      {
	 vector_inner_prod(M, A+i, 1, B+i, 1, Acc[i]);
      }
      vector_sum(N, Acc.cptr(), 1, r);
   }
   else
   {
      PANIC("not implemented");
   }
}

void
matrix_inner_prod(char Atrans, char Btrans, int M, int N,
		  cuda::const_gpu_ptr<std::complex<double>> A, int ldA,
		  cuda::const_gpu_ptr<std::complex<double>> B, int ldB,
		  cuda::gpu_ref<std::complex<double>>& r)
{
   if (Atrans == 'N' && Btrans == 'N')
   {
      cuda::gpu_buffer<std::complex<double>> Acc = cuda::allocate_gpu_temporary<std::complex<double>>(N);

      for (int i = 0; i < N; ++i)
      {
	 vector_inner_prod(M, A+i, 1, B+i, 1, Acc[i]);
      }
      vector_sum(N, Acc.cptr(), 1, r);
   }
   else
   {
      PANIC("not implemented");
   }
}

void
matrix_add_inner_prod(char Atrans, char Btrans, int M, int N,
		  cuda::const_gpu_ptr<std::complex<double>> A, int ldA,
		  cuda::const_gpu_ptr<std::complex<double>> B, int ldB,
		  cuda::gpu_ref<std::complex<double>>& r)
{
   if (Atrans == 'N' && Btrans == 'N')
   {
      gpu_buffer<std::complex<double>> Acc = cuda::allocate_gpu_temporary<std::complex<double>>(N+1);
      Acc[0] = r;

      for (int i = 0; i < N; ++i)
      {
	 vector_inner_prod(M, A+i, 1, B+i, 1, Acc[i]);
      }
      vector_sum(N, Acc.cptr(), 1, r);
   }
   else
   {
      PANIC("not implemented");
   }
}

} // namespace cuda
