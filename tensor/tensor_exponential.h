// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensor_exponential.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_TENSOR_TENSOR_EXPONENTIAL_H)
#define MPTOOLKIT_TENSOR_TENSOR_EXPONENTIAL_H

#include "tensor.h"
#include "blas/matrix.h"

namespace Tensor
{

// exponentiate a scalar operator
inline
IrredTensor<std::complex<double>, BasisList, BasisList>
exp(IrredTensor<std::complex<double>, BasisList, BasisList> const& m);

IrredTensor<blas::Matrix<std::complex<double>, blas::cpu_tag>, VectorBasis, VectorBasis>
exp(IrredTensor<blas::Matrix<std::complex<double>, blas::cpu_tag>, VectorBasis, VectorBasis> const& m);

// the exponential is done on the CPU, so with a device tag we simply copy to/from the CPU
template <typename Tag>
IrredTensor<blas::Matrix<std::complex<double>, Tag>, VectorBasis, VectorBasis>
exp(IrredTensor<blas::Matrix<std::complex<double>, Tag>, VectorBasis, VectorBasis> const& m);

} // namespace Tensor

#endif
