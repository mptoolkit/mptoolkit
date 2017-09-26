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

// intializes cublas to run in a thread - must be called once per thread prior to
// making any other cublas calls.  The CUDA device must be intialized prior to this call.
void setup_cublas_thread();

// returns the thread-local handle
handle& get_handle();

} // namespace cublas

#include "cublas.icc"

#endif
