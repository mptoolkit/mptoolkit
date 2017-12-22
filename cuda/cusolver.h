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

char const* GetErrorString(cusolverStatus_t error);

// wrapper for a cusolverStatus_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cusolverStatus_t Err) : std::runtime_error(std::string("cuSOLVER error: ")
                                                       + cusolver::GetErrorString(Err)), err_(Err)
      {std::cerr << "cuSOLVER Error " << int(Err) << ' ' << cusolver::GetErrorString(Err) << '\n';}

      cusolverStatus_t code() const { return err_; }
      operator cusolverStatus_t() const { return err_; }

      char const* name() const { return cusolver::GetErrorName(err_); }
      char const* string() const { return cusolver::GetErrorString(err_); }

   private:
      cusolverStatus_t err_;
};

inline
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

      static handle create() { cusolverDnHandle_t h; cusolverDnCreate(&h); return handle(h); }

      void destroy() { cusolverDnDestroy(h_); h_ = nullptr; }

      // set the stream associated with the handle
      void set_stream(cuda::stream const& s)
      {
         check_error(cusolverDnSetStream(h_, s.raw_stream()));
      }

   private:
      handle(cusolverDnHandle_t h) : h_(h) {}

      cusolverDnHandle_t h_;
};

// returns the thread-local handle
handle& get_handle();

} // namespace cusolver


namespace cuda
{

// LAPACK functions must go in namespace cuda so they are found during ADL

void DiagonalizeSymmetric(int Size, cuda::const_gpu_ptr<double> A, int ldA, cuda::gpu_ptr<double> Eigen);

void DiagonalizeHermitian(int Size, cuda::const_gpu_ptr<std::complex<double>> A, int ldA, 
			  cuda::gpu_ptr<double> Eigen);

void SingularValueDecomposition(int Rows, int Cols, 
				cuda::gpu_ptr<double> Data, int LeadingDim, 
				cuda::gpu_ptr<double> Dvec,
				cuda::gpu_ptr<double> Umat, int ldU, 
				cuda::gpu_ptr<double> Vmat, int ldV);

void SingularValueDecomposition(int Rows, int Cols, 
				cuda::gpu_ptr<std::complex<double>> Data, int LeadingDim, 
				cuda::gpu_ptr<double> Dvec,
				cuda::gpu_ptr<std::complex<double>> Umat, int ldU, 
				cuda::gpu_ptr<std::complex<double>> Vmat, int ldV);

} // namespace cuda

#include "cusolver.icc"

#endif
