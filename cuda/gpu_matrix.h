// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cublas/gpu_matrix.h
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
// gpu_matrix is a matrix stored in device memory on the GPU.
// The basic operations are implemented by cuBLAS calls, which are
// non-blocking.  Blocking calls have a wait_ prefix.
//
// LAPACK functionality (eg implemented by cuSOLVER) are generally
// blocking calls.  There is no attempt to make these asynchronous.
// Asynchronous LAPACK calls will be handled by a higher level library.
//

#if !defined(MPTOOLKIT_CUBLAS_GPU_MATRIX_H)
#define MPTOOLKIT_CUBLAS_GPU_MATRIX_H

#include "blas/matrixref.h"
#include "blas/matrix.h"
#include "blas/vector_view.h"
#include "cublas.h"
#include "cusolver.h"
#include "gpu_vector.h"
#include "gpu_buffer.h"

namespace blas
{

template <typename T>
using gpu_matrix = blas::Matrix<T, gpu_tag>;

// blocking matrix get
template <typename T>
blas::Matrix<T>
get_wait(gpu_matrix<T> const& M)
{
   blas::Matrix<T> Result(M.rows(), M.cols());
   cublas::check_error(cublasGetMatrix(M.rows(), M.cols(), sizeof(T),
				       M.storage().device_ptr(), M.leading_dimension(),
				       Result.storage(), Result.leading_dimension()));
   return Result;
}

// blocking matrix set
template <typename T>
void
set_wait(gpu_matrix<T>& A, blas::Matrix<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   //TRACE(A.cols())(A.rows())(A.device_ptr())(A.leading_dim())(B.data())(leading_dimension(B));
   cublas::check_error(cublasSetMatrix(A.rows(), A.cols(), sizeof(T),
				       B.data(), leading_dimension(B),
				       A.device_ptr(), A.leading_dim()));
}

// non-blocking set
template <typename T>
cuda::event
set(gpu_matrix<T>& A, blas::Matrix<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   cublas::setMatrixAsync(A.rows(), A.cols(), B.storage(), B.leading_dimension(),
                          A.storage(), A.leading_dimension());
   return A.storage().sync();
}

} // namespace cublas

#endif
