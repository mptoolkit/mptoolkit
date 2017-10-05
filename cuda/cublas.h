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

#if !defined(MPTOOLKIT_CUDA_CUBLAS_H)
#define MPTOOLKIT_CUDA_CUDBLAS_H

#include "cuda.h"
#include "gpu_buffer.h"
#include <list>
#include <mutex>
#include <cublas_v2.h>

#include <iostream>

namespace cublas
{

// returns the cublas version number
int version();

// TODO: cublas error class

inline
char const* cublasGetErrorName(cublasStatus_t error)
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

inline
char const* cublasGetErrorString(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "cuBLAS library not initialized";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "memory allocation failed";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "invalid value or parameter";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "feature not supported by this architecture (possible double-precision?)";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "invalid GPU memory mapping";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "GPU kernel execution failed";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "internal error";
    }

    return "<cublas-unknown>";
}

inline
void check_error(cublasStatus_t s)
{
   if (s != CUBLAS_STATUS_SUCCESS)
   {
      std::string ss = std::string("cuBLAS error: ") + cublasGetErrorName(s);
      throw std::runtime_error(ss);
   }
}

// cublas handle.  Moveable, but not copyable.
class handle
{
   public:
      handle() : h_(nullptr) {}
      handle(handle&& other) : h_(other.h_) { other.h_ = nullptr; }
      handle(handle const&) = delete;
      handle& operator=(handle&& other) { std::swap(h_, other.h_); return *this; }
      handle& operator=(handle const&) = delete;
      ~handle() { if (h_) cublasDestroy(h_); }

      cublasHandle_t raw_handle() const { return h_; }

      static handle create() { cublasHandle_t h; cublasCreate(&h); return handle(h); }

      void destroy() { cublasDestroy(h_); h_ = nullptr; }

      // set the stream associated with the handle
      void set_stream(cuda::stream const& s)
      {
         check_error(cublasSetStream(h_, s.raw_stream()));
      }

      void set_pointer_mode(cublasPointerMode_t m)
      {
         check_error(cublasSetPointerMode(h_, m));
      }

   private:
      handle(cublasHandle_t h) : h_(h) {}

      cublasHandle_t h_;
};

inline
cublasOperation_t cublas_trans(char c)
{
   switch (c)
   {
   case 'N' : return CUBLAS_OP_N;
   case 'T' : return CUBLAS_OP_T;
   case 'C' : return CUBLAS_OP_C;
   default : throw std::runtime_error("Unsupported TRANS in cublas");
   }
}

// intializes cublas to run in a thread - must be called once per thread prior to
// making any other cublas calls.  The CUDA device must be intialized prior to this call.
void setup_cublas_thread();

// returns the thread-local handle
handle& get_handle();

inline
void
gemm(cublas::handle& H, char Atrans, char Btrans, int M, int N, int K, double alpha,
     cuda::gpu_buffer<double> const& A, int lda, cuda::gpu_buffer<double> const& B, int ldb,
     double beta, cuda::gpu_buffer<double>& C, int ldc)
{
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   check_error(cublasDgemm(H.raw_handle(), cublas_trans(Atrans), cublas_trans(Btrans), M, N, K,
                           &alpha, A.device_ptr(), lda, B.device_ptr(), ldb,
                           &beta, C.device_ptr(), ldc));
}

// geam - we have two versions, for in-place and out-of-place operations
inline
void
geam(cublas::handle& H, char Atrans, int M, int N, double alpha,
     cuda::gpu_buffer<double> const& A, int lda, cuda::gpu_buffer<double>& C, int ldc)
{
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   double beta = 0.0;
   check_error(cublasDgeam(H.raw_handle(), cublas_trans(Atrans), CUBLAS_OP_N, M, N,
                           &alpha, A.device_ptr(), lda,
                           &beta, C.device_ptr(), ldc,
                           C.device_ptr(), ldc));
}

inline
void
geam(cublas::handle& H, char Atrans, char Btrans, int M, int N,
     double alpha, cuda::gpu_buffer<double> const& A, int lda,
     double beta,  cuda::gpu_buffer<double> const& B, int ldb,
     cuda::gpu_buffer<double>& C, int ldc)
{
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   check_error(cublasDgeam(H.raw_handle(), cublas_trans(Atrans), cublas_trans(Btrans), M, N,
                           &alpha, A.device_ptr(), lda,
                           &beta, B.device_ptr(), ldb,
                           C.device_ptr(), ldc));
}

template <typename T>
inline
void
setMatrixAsync(int N, int M, T const* A, int lda, cuda::gpu_buffer<T>& B, int ldb)
{
   check_error(cublasSetMatrixAsync(N, M, sizeof(T), A, lda, B.device_ptr(), ldb,
                                    B.get_stream().raw_stream()));
}

} // namespace cublas

#include "cublas.icc"

#endif
