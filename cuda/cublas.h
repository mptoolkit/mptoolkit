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
#define MPTOOLKIT_CUDA_CUBLAS_H

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

// helper functions to get a name and error description from a cublasStatus_t
char const* cublasGetErrorName(cublasStatus_t error);

char const* cublasGetErrorString(cublasStatus_t error);

// wrapper for a cublasStatus_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cublasStatus_t Err) : std::runtime_error(std::string("cuBLAS error: ") + cublasGetErrorString(Err)), err_(Err)
      {std::cerr << "cuBLAS Error " << int(Err) << ' ' << cublasGetErrorString(Err) << '\n';}

      cublasStatus_t code() const { return err_; }
      operator cublasStatus_t() const { return err_; }

      char const* name() const { return cublasGetErrorName(err_); }
      char const* string() const { return cublasGetErrorString(err_); }

   private:
      cublasStatus_t err_;
};

inline
void check_error(cublasStatus_t s)
{
   if (s != CUBLAS_STATUS_SUCCESS)
   {
      throw error(s);
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

// utility functions

template <typename T>
inline
void
setMatrixAsync(int N, int M, T const* A, int lda, cuda::gpu_ptr<T> B, int ldb)
{
   check_error(cublasSetMatrixAsync(N, M, sizeof(T), A, lda, B.device_ptr(), ldb,
                                    B.get_stream().raw_stream()));
}

template <typename T>
inline
void
setVectorAsync(int N, T const* A, int stridea, cuda::gpu_ptr<T> B, int strideb)
{
   check_error(cublasSetVectorAsync(N, sizeof(T), A, stridea, B.device_ptr(), strideb,
                                    B.get_stream().raw_stream()));
}

// BLAS level 1

inline
void
copy(cublas::handle& H, int n, cuda::const_gpu_ptr<double> x, int incx,
     cuda::gpu_ptr<double> y, int incy)
{
   H.set_stream(y.get_stream());
   y.wait_for(x);
   check_error(cublasDcopy(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy));
}

inline
void
scal(cublas::handle& H, int n, double alpha, cuda::gpu_ptr<double> y, int incy)
{
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   check_error(cublasDscal(H.raw_handle(), n, &alpha, y.device_ptr(), incy));
}

inline
void
axpy(cublas::handle& H, int n, double alpha,
     cuda::const_gpu_ptr<double> x, int incx,
     cuda::gpu_ptr<double> y, int incy)
{
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   y.wait_for(x);
   check_error(cublasDaxpy(H.raw_handle(), n, &alpha, x.device_ptr(), incx, y.device_ptr(), incy));
}

// BLAS level 2

inline
void
gemv(cublas::handle& H, char Atrans, int M, int N, double alpha,
     cuda::const_gpu_ptr<double> A, int lda, cuda::const_gpu_ptr<double> x, int incx,
     double beta, cuda::gpu_ptr<double> y, int incy)
{
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   y.wait_for(A);
   y.wait_for(x);
   check_error(cublasDgemv(H.raw_handle(), cublas_trans(Atrans), M, N,
                           &alpha, A.device_ptr(), lda, x.device_ptr(), incx,
                           &beta, y.device_ptr(), incy));
}

// BLAS level 2 extensions

// geam - we have two versions, for in-place and out-of-place operations
inline
void
geam(cublas::handle& H, char Atrans, int M, int N, double alpha,
     cuda::const_gpu_ptr<double> A, int lda, cuda::gpu_ptr<double> C, int ldc)
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
     double alpha, cuda::const_gpu_ptr<double> A, int lda,
     double beta,  cuda::const_gpu_ptr<double> B, int ldb,
     cuda::gpu_ptr<double> C, int ldc)
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

// BLAS level 3

inline
void
gemm(cublas::handle& H, char Atrans, char Btrans, int M, int N, int K, double alpha,
     cuda::const_gpu_ptr<double> A, int lda, cuda::const_gpu_ptr<double> B, int ldb,
     double beta, cuda::gpu_ptr<double> C, int ldc)
{
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   check_error(cublasDgemm(H.raw_handle(), cublas_trans(Atrans), cublas_trans(Btrans), M, N, K,
                           &alpha, A.device_ptr(), lda, B.device_ptr(), ldb,
                           &beta, C.device_ptr(), ldc));
}

} // namespace cublas

#include "cublas.icc"

#endif
