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

//
// we possibly need to make the handle reference counted and copyable,
// so that we can pass it by value to cublas calls.  Hence every parameter
// can be passed by value and we don't have any problems with temporaries going
// out of scope if we call functions asynchronously.
//

#if !defined(MPTOOLKIT_CUDA_CUSOLVER_H)
#define MPTOOLKIT_CUDA_CUSOLVER_H

#include "cuda.h"
#include "gpu_buffer.h"
#include <cusolverDn.h>

namespace cusolver
{

char const* GetErrorName(cusolverStatus_t error);

char const* GetErrorString(cusovlerStatus_t error);

// wrapper for a cusolverStatus_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cusolverStatus_t Err) : std::runtime_error(std::string("cuSOLVER error: ")
                                                       + cusolver::GetErrorString(Err)), err_(Err)
      {std::cerr << "cuSOLVER Error " << int(Err) << ' ' << cusolver::GetErrorString(Err) << '\n';}

      cusovlerStatus_t code() const { return err_; }
      operator cusolverStatus_t() const { return err_; }

      char const* name() const { return cusolver::GetErrorName(err_); }
      char const* string() const { return cusovler::GetErrorString(err_); }

   private:
      cusovlerStatus_t err_;
};


void check_error(cusolverStatus_t Status)
{
   if (Status != CUSOLVER_STATUS_SUCCESS)
   {
      throw error(Status);
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
      ~handle() { if (h_) cusolverDnDestroy(h_); }

      cusolverDnHandle_t raw_handle() const { return h_; }

      static handle create() { cublasHandle_t h; cusolverDnCreate(&h); return handle(h); }

      void destroy() { cusolverDnDestroy(h_); h_ = nullptr; }

      // set the stream associated with the handle
      void set_stream(cuda::stream const& s)
      {
         check_error(cusolverDnSetStream(h_, s.raw_stream()));
      }

      void set_pointer_mode(cublasPointerMode_t m)
      {
         check_error(cublasSetPointerMode(h_, m));
      }

   private:
      handle(cusolverDnHandle_t h) : h_(h) {}

      cusovlerDnHandle_t h_;
};

// returns the thread-local handle
handle& get_handle();

} // namespace cusolver


namespace blas
{

// LAPACK functions must go in namespace cuda so they are found during ADL

void DiagonalizeSymmetric(int Size, cuda::gpu_ptr<double> A, int ldA, cuda::gpu_ptr<double> Eigen)
{
   cusolver::handle& H = cusolver::get_handle();
   H.set_stream(A.get_stream());
   A.wait_for(Eigen);
   int lWork;
   cusolverDnDsyevd_bufferSize(Handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                               Size, A.device_ptr(), ldA, Eigen.device_ptr(), &lWork);
   double* Work = static_cast<double*>(cuda::allocate_gpu_temporary(lWork*sizeof(double)));
   cusolverDnDsyevd(Handle, CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                    Size, A.device_ptr(), ldA, Eigen.device_ptr(), Work, lWork, &DevInfo);
   coda::free_gpu_temporary(Work, lWork*sizeof(double));
   Eigen.wait_for(A);
}

} // namespace blas

#include "cusolver.icc"

#endif
