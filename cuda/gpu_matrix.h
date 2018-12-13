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
#include "blas/diagonalmatrix.h"
#include "cublas.h"
#include "cusolver.h"
#include "gpu_vector.h"
#include "gpu_buffer.h"

namespace blas
{

template <typename T>
using gpu_matrix = Matrix<T, gpu_tag>;

template <typename T>
using gpu_diagonal_matrix = DiagonalMatrix<T, gpu_tag>;


// blocking matrix get
template <typename T>
Matrix<T>
get_wait(gpu_matrix<T> const& M)
{
   M.storage().get_stream().synchronize(); // Wait until pending operations on M are complete
   Matrix<T> Result(M.rows(), M.cols());
   cublas::check_error(cublasGetMatrix(M.rows(), M.cols(), sizeof(T),
				       M.storage().device_ptr(), M.leading_dimension(),
				       Result.storage(), Result.leading_dimension()));
   return Result;
}

// blocking matrix set
template <typename T>
void
set_wait(gpu_matrix<T>& A, Matrix<T> const& B)
{
   A.storage().get_stream().synchronize(); // Wait until pending operations on A are complete, otherwise they will overwrite the matrix
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   //TRACE(A.cols())(A.rows())(A.device_ptr())(A.leading_dim())(B.data())(leading_dimension(B));
   cublas::check_error(cublasSetMatrix(A.rows(), A.cols(), sizeof(T),
				       B.storage(), B.leading_dimension(),
				       A.storage().device_ptr(), A.leading_dimension()));
}

// non-blocking set
template <typename T>
cuda::event
set(gpu_matrix<T>& A, Matrix<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   cublas::setMatrixAsync(A.rows(), A.cols(), B.storage(), B.leading_dimension(),
                          A.storage(), A.leading_dimension());
   return A.storage().sync();
}

template <typename T>
std::ostream& operator<<(std::ostream& out, gpu_matrix<T> const& A)
{
   out << "gpu_matrix<" << tracer::typeid_name<T>() << "> ";
   //#if !defined(NDEBUG)
   out << get_wait(A);
   //#else
   // out << '[' << A.rows() << "," << A.cols() << ']';
   //#endif
   return out;
}

// DiagonalMatrix

template <typename T>
DiagonalMatrix<T>
get_wait(gpu_diagonal_matrix<T> const& M)
{
   DiagonalMatrix<T> Result(M.rows());
   Result.diagonal() = get_wait(M.diagonal());
   return Result;
}

template <typename T>
void
set_wait(gpu_diagonal_matrix<T>& A, DiagonalMatrix<T> const& B)
{
   DEBUG_CHECK_EQUAL(A.rows(), B.rows());
   DEBUG_CHECK_EQUAL(A.cols(), B.cols());
   set_wait(A.diagonal(), B.diagonal());
}

// TODO: non-blocking

} // namespace cublas

#endif
