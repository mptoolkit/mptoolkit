// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/detail/cublas_detail.h
//
// Copyright (C) 2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_CUDA_DETAIL_CUSOLVER_DETAIL_H)
#define MPTOOLKIT_CUDA_DETAIL_CUSOLVER_DETAIL_H

#include "cuda/cuda.h"
#include "cuda/gpu_buffer.h"
#include "cuda/gpu_vector.h"

namespace blas
{

namespace detail
{

// LAPACK functions must go in namespace blas::detail

void DiagonalizeSymmetric(int Size, cuda::gpu_ptr<double> A, int ldA, cuda::gpu_ptr<double> Eigen);

void DiagonalizeHermitian(int Size, cuda::gpu_ptr<std::complex<double>> A, int ldA,
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

// for MxN matrix, if N>M then fill all rows of Vmat.  If M<N then fill all columns of Umat.
// Note that the extra zero entries to pad out Dvec and UMmat (respectively Vmat) are NOT
// set by this function; these need to be set by the caller.
void SingularValueDecompositionFull(int Rows, int Cols,
				    cuda::gpu_ptr<double> Data, int LeadingDim,
				    cuda::gpu_ptr<double> Dvec,
				    cuda::gpu_ptr<double> Umat, int ldU,
				    cuda::gpu_ptr<double> Vmat, int ldV);

void SingularValueDecompositionFull(int Rows, int Cols,
				    cuda::gpu_ptr<std::complex<double>> Data, int LeadingDim,
				    cuda::gpu_ptr<double> Dvec,
				    cuda::gpu_ptr<std::complex<double>> Umat, int ldU,
				    cuda::gpu_ptr<std::complex<double>> Vmat, int ldV);

} // namespace detail

} // namespace blas

#endif
