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
#include "gpu_vector.h"
#include "detail/cusolver_detail.h"
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
      throw error(Status);
#if defined(CUDA_SYNCHRONIZE)
   cudaError_t e = cudaDeviceSynchronize();
   if (e != cudaSuccess)
      throw cuda::error(e);
#endif
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

// ARGH, CUBLAS svd requires transpositions in some cases,
// so overload the blas SVD functions for gpu buffers.
// This isn't a great solution - better would be to implement
// the transpositions in detail/cusolver_detail.h, but that is more
// complicated since those functions only have access to data buffers,
// not actual matrices.

namespace blas
{

struct gpu_tag;

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V>
void
SVD(NormalMatrix<Scalar, M, gpu_tag>&& Mmat, NormalMatrix<Scalar, U, blas::gpu_tag>& Umat,
    NormalVector<Real, D, gpu_tag>& Dvec, NormalMatrix<Scalar, V, blas::gpu_tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), std::min(Mmat.rows(), Mmat.cols()));
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());

   if (Mmat.rows() < Mmat.cols())
   {
      Matrix<Scalar, gpu_tag> TempM(trans(Mmat));
      Matrix<Scalar, gpu_tag> TempU(Vmat.cols(), Vmat.rows());
      Matrix<Scalar, gpu_tag> TempV(Umat.cols(), Umat.rows());
      SVD(std::move(TempM), TempU, Dvec, TempV);
      assign(Umat, trans(std::move(TempV)));
      assign(Vmat, trans(std::move(TempU)));
   }
   else
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
                                         Umat.storage(), Umat.leading_dimension(),
                                         Vmat.storage(), Vmat.leading_dimension());
   }
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V>
void
SVD(NormalMatrix<Scalar, M, gpu_tag>&& Mmat, NormalMatrix<Scalar, U, gpu_tag>& Umat,
    NormalVectorProxy<Real, D, gpu_tag>&& Dvec, NormalMatrix<Scalar, V, gpu_tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), std::min(Mmat.rows(), Mmat.cols()));
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());

   if (Mmat.rows() < Mmat.cols())
   {
      Matrix<Scalar, gpu_tag> TempM(trans(Mmat));
      Matrix<Scalar, gpu_tag> TempU(Vmat.cols(), Vmat.rows());
      Matrix<Scalar, gpu_tag> TempV(Umat.cols(), Umat.rows());
      SVD(std::move(TempM), TempU, std::move(Dvec), TempV);
      assign(Umat, trans(std::move(TempV)));
      assign(Vmat, trans(std::move(TempU)));
   }
   else
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
                                         std::move(Dvec).storage(),
                                         Umat.storage(), Umat.leading_dimension(),
                                         Vmat.storage(), Vmat.leading_dimension());
   }
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V>
void
SVD_FullRows(NormalMatrix<Scalar, M, gpu_tag>&& Mmat, NormalMatrix<Scalar, U, gpu_tag>& Umat,
	     NormalVector<Real, D, gpu_tag>& Dvec, NormalMatrix<Scalar, V, gpu_tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.rows());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.rows() < Mmat.cols())
   {
      Matrix<Scalar, gpu_tag> TempM(trans(Mmat));
      Matrix<Scalar, gpu_tag> TempU(Vmat.cols(), Vmat.rows());
      Matrix<Scalar, gpu_tag> TempV(Umat.cols(), Umat.rows());
      SVD_FullCols(std::move(TempM), TempU, Dvec, TempV);
      assign(Umat, trans(std::move(TempV)));
      assign(Vmat, trans(std::move(TempU)));
   }
   else
   {
      // need to zero the additional elements in Dvec, for remaining zero singular values beyond the number of rows
      // clear(Dvec);
      // not needed now, SingularValueDecompositionFull zeroes out the elements
      detail::SingularValueDecompositionFull(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
					     Umat.storage(), Umat.leading_dimension(),
					     Vmat.storage(), Vmat.leading_dimension());
   }
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V>
void
SVD_FullRows(NormalMatrix<Scalar, M, gpu_tag>&& Mmat, NormalMatrix<Scalar, U, gpu_tag>& Umat,
	     NormalVectorProxy<Real, D, gpu_tag>&& Dvec, NormalMatrix<Scalar, V, gpu_tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.rows());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.rows() < Mmat.cols())
   {
      Matrix<Scalar, gpu_tag> TempM(trans(Mmat));
      Matrix<Scalar, gpu_tag> TempU(Vmat.cols(), Vmat.rows());
      Matrix<Scalar, gpu_tag> TempV(Umat.cols(), Umat.rows());
      SVD_FullCols(std::move(TempM), TempU, std::move(Dvec), TempV);
      assign(Umat, trans(std::move(TempV)));
      assign(Vmat, trans(std::move(TempU)));
   }
   else
   {
      // need to zero the additional elements in Dvec
      // not needed now, SingularValueDecompositionFull zeroes out the elements
      // clear(Dvec);
      detail::SingularValueDecompositionFull(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
					     std::move(Dvec).storage(),
					     Umat.storage(), Umat.leading_dimension(),
					     Vmat.storage(), Vmat.leading_dimension());
      DEBUG_CHECK(!isnan(norm_frob(Umat)))(Umat)(Mmat);
      DEBUG_CHECK(!isnan(norm_frob(Vmat)))(Vmat)(Mmat);
   }
}

// SVD_FullCols
//
// Given M as an m*n matrix, U is m*n, D is n*n, V is n*n
//

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V>
void
SVD_FullCols(NormalMatrix<Scalar, M, gpu_tag>&& Mmat, NormalMatrix<Scalar, U, gpu_tag>& Umat,
	     NormalVector<Real, D, gpu_tag>& Dvec, NormalMatrix<Scalar, V, gpu_tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.cols());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.cols() <= Mmat.rows())
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(),
                                         Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
					 Umat.storage(), Umat.leading_dimension(),
					 Vmat.storage(), Vmat.leading_dimension());
   }
   else
   {
      Matrix<Scalar, gpu_tag> TempM(trans(Mmat));
      Matrix<Scalar, gpu_tag> TempU(Vmat.cols(), Vmat.rows());
      Matrix<Scalar, gpu_tag> TempV(Umat.cols(), Umat.rows());
      SVD_FullRows(std::move(TempM), TempU, Dvec, TempV);
      assign(Umat, trans(std::move(TempV)));
      assign(Vmat, trans(std::move(TempU)));
   }
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V>
void
SVD_FullCols(NormalMatrix<Scalar, M, blas::gpu_tag>&& Mmat, NormalMatrix<Scalar, U, gpu_tag>& Umat,
	     NormalVectorProxy<Real, D, blas::gpu_tag>&& Dvec, NormalMatrix<Scalar, V, gpu_tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.cols());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.cols() <= Mmat.rows())
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
					 std::move(Dvec).storage(),
					 Umat.storage(), Umat.leading_dimension(),
					 Vmat.storage(), Vmat.leading_dimension());
   }
   else
   {
      Matrix<Scalar, gpu_tag> TempM(trans(Mmat));
      Matrix<Scalar, gpu_tag> TempU(Vmat.cols(), Vmat.rows());
      Matrix<Scalar, gpu_tag> TempV(Umat.cols(), Umat.rows());
      SVD_FullRows(std::move(TempM), TempU, std::move(Dvec), TempV); 
      DEBUG_CHECK(!isnan(norm_frob(TempU)))(TempU);
      DEBUG_CHECK(!isnan(norm_frob(TempV)))(TempV);
      assign(Umat, trans(std::move(TempV)));
      assign(Vmat, trans(std::move(TempU)));
   }
}

} // namespace blas

#include "cusolver.icc"

#endif
