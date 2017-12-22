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
#include "blas/diagonalmatrix.h"

using namespace Tensor;

// matrix types

namespace cpu
{
using Matrix             = blas::Matrix<complex>;
using RealMatrix         = blas::Matrix<real>;
using DiagonalMatrix     = blas::DiagonalMatrix<complex>;
using RealDiagonalMatrix = blas::DiagonalMatrix<real>;
using Vector             = blas::Vector<complex>;
using RealVector         = blas::Vector<real>;
} // namspace cpu

#if defined(HAVE_CUDA)
#include "cuda/gpu_matrix.h"
#include "cuda/gpu_vector.h"

using device_tag = blas::gpu_tag;

#elif defined(USE_THREADS)
#error "thread pool not implemented yet"

#else
// Fallback to synchronous CPU

using device_tag = blas::cpu_tag;

#endif

using Matrix             = blas::Matrix<complex, device_tag>;
using RealMatrix         = blas::Matrix<real, device_tag>;
using DiagonalMatrix     = blas::DiagonalMatrix<complex, device_tag>;
using RealDiagonalMatrix = blas::DiagonalMatrix<real, device_tag>;
using Vector             = blas::Vector<complex, device_tag>;
using RealVector         = blas::Vector<real, device_tag>;

// Tensor types.  The matrix versions use device matrices

using SimpleOperator       = IrredTensor<complex>;
using RealSimpleOperator   = IrredTensor<real>;

using MatrixOperator       = IrredTensor<Matrix,
                                         VectorBasis,
                                         VectorBasis>;

using RealMatrixOperator   = IrredTensor<RealMatrix,
                                         VectorBasis,
                                         VectorBasis>;

using RealDiagonalOperator = IrredTensor<RealDiagonalMatrix,
                                         VectorBasis,
                                         VectorBasis>;//,
//                                         Tensor::DiagonalStructure>;

using SimpleRedOperator = ReducibleTensor<complex, BasisList, BasisList>;

#endif
