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
IrredTensor<std::complex<double>, BasisList, BasisList>
exp(IrredTensor<std::complex<double>, BasisList, BasisList> const& m);

template <typename T, typename Tag>
IrredTensor<blas::Matrix<T, Tag>, VectorBasis, VectorBasis>
exp(IrredTensor<blas::Matrix<T, Tag>, VectorBasis, VectorBasis> const& m);

} // namespace Tensor

#include "tensor_exponential.icc"

#endif
