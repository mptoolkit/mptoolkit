// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/types.h
//
// Copyright (C) 2019 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_TENSOR_TYPES_H)
#define MPTOOLKIT_TENSOR_TYPES_H

#include "common/types.h"
#include "tensor/tensor.h"
#include "tensor/reducible.h"
#include "blas/matrix.h"
#include "blas/matrix-eigen.h"

using namespace Tensor;

// matrix types

using Matrix             = blas::Matrix<complex>;
using RealMatrix         = blas::Matrix<real>;
using DiagonalMatrix     = blas::DiagonalMatrix<complex>;
using RealDiagonalMatrix = blas::DiagonalMatrix<real>;
using Vector             = blas::Vector<complex>;
using RealVector         = blas::Vector<real>;

#if defined(HAVE_CUDA)
#include "cuda/gpu_matrix.h"
#include "cuda/gpu_vector.h"
//#include "cuda/gpu_eigen.h"

using Matrix_Device             = cuda::gpu_matrix<complex>;
using RealMatrix_Device         = cuda::gpu_matrix<real>;
using DiagonalMatrix_Device     = cuda::gpu_diagonal_matrix<complex>;
using RealDiagonalMatrix_Device = cuda::gpu_diagonal_matrix<real>;
using Vector_Device             = cuda::gpu_vector<complex>;
using RealVector_Device         = cuda::gpu_vector<real>;

#elif defined(USE_THREADS)
#error "thread pool not implemented yet"

#else
// Fallback to synchronous CPU

using Matrix_Device = Matrix;
using RealMatrix_Device = RealMatrix;
using DiagonalMatrix_Device = DiagonalMatrix;
using RealDiagonalMatrix_Device = RealDiagonalMatrix;
using Vector_Device = Vector;
using RealVector_Device = RealVector;

#endif

// Tensor types.  The matrix versions use device matrices

using SimpleOperator       = IrredTensor<complex>;
using RealSimpleOperator   = IrredTensor<real>;

using MatrixOperator       = IrredTensor<LinearAlgebra::Matrix_Device<complex>,
                                         VectorBasis,
                                         VectorBasis>
using RealMatrixOperator   = IrredTensor<LinearAlgebra::Matrix_Device<real>,
                                         VectorBasis,
                                         VectorBasis>

using RealDiagonalOperator = IrredTensor<LinearAlgebra::DiagonalMatrix_Device<double>,
                                         VectorBasis,
                                         VectorBasis,
                                         Tensor::DiagonalStructure>;

#endif
