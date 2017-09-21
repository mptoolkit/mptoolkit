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
#include "linearalgebra2/matrix.h"
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

   private:
      handle(cublasHandle_t h) : h_(h) {}

      cublasHandle_t h_;
};

// set the stream associated with the given cublas handle
void set_stream(handle const& h, cuda::stream const& s);

inline
void set_stream(handle const& h, cuda::stream const& s)
{
   check_error(cublasSetStream(h.raw_handle(), s.raw_stream()));
}

//
// GPU-storage matrix type.  Column-major format for compatability with BLAS
//

template <typename T>
class gpu_matrix
{
   public:
      typedef T value_type;

      gpu_matrix(int Rows_, int Cols_, cuda::arena const& A);

      int rows() const { return Rows; }
      int cols() const { return Cols; }

      int leading_dim() const { return LeadingDimension; }

      T* device_ptr() { return Buf.device_ptr(); }
      T const* device_ptr() const { return Buf.device_ptr(); }

   private:
      int Rows;
      int Cols;
      int LeadingDimension;
      cuda::gpu_buffer<T> Buf;
};

template <typename T>
inline
gpu_matrix<T>::gpu_matrix(int Rows_, int Cols_, cuda::arena const& A)
   : Rows(Rows_), Cols(Cols_), LeadingDimension(Rows_),
     Buf(cuda::gpu_buffer<T>::allocate(Rows_*Cols_, A))
{
}

// blocking matrix get
template <typename T>
LinearAlgebra::Matrix<T>
get_wait(gpu_matrix<T> const& M)
{
   LinearAlgebra::Matrix<T> Result(M.rows(), M.cols());
   cublas::check_error(cublasGetMatrix(M.rows(), M.cols(), sizeof(T),
				       M.device_ptr(), M.leading_dim(),
				       Result.data(), leading_dimension(Result)));
   return Result;
}

template <typename T>
void
set_wait(gpu_matrix<T>& A, LinearAlgebra::Matrix<T> const& B)
{
   //TRACE(A.cols())(A.rows())(A.device_ptr())(A.leading_dim())(B.data())(leading_dimension(B));
   cublas::check_error(cublasSetMatrix(A.rows(), A.cols(), sizeof(T),
				       B.data(), leading_dimension(B),
				       A.device_ptr(), A.leading_dim()));
}

// generic gemm = C' = alpha*A*B + beta*C
inline
void gemm(handle& H, double alpha, gpu_matrix<double> const& A, gpu_matrix<double> const& B, double beta, gpu_matrix<double>& C)
{
   check_error(cublasDgemm(H.raw_handle(), CUBLAS_OP_N, CUBLAS_OP_N, A.rows(), A.cols(), B.cols(), 
			    &alpha, A.device_ptr(), A.leading_dim(), B.device_ptr(), B.leading_dim(),
			    &beta, C.device_ptr(), C.leading_dim()));
}

} // namespace cublas

#endif
