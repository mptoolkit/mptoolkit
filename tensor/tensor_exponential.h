// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/tensor_exponential.h
//
// Copyright (C) 2004-2020 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(TENSOR_EXPONENTIAL_H_DHCKJDHUREYT7845Y78Y78TY78TY78T)
#define TENSOR_EXPONENTIAL_H_DHCKJDHUREYT7845Y78Y78TY78TY78T

#include "tensor.h"

namespace Tensor
{

// exponentiate a scalar operator
IrredTensor<std::complex<double>, BasisList, BasisList>
Exponentiate(IrredTensor<std::complex<double>, BasisList, BasisList> const& m);

inline
IrredTensor<std::complex<double>, BasisList, BasisList>
exp(IrredTensor<std::complex<double>, BasisList, BasisList> const& m)
{
   return Exponentiate(m);
}

IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
exp(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> const& m);

IrredTensor<std::complex<double>, BasisList, BasisList>
abs(IrredTensor<std::complex<double>, BasisList, BasisList> const& m);

IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
abs(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> const& m);

} // namespace Tensor

#endif
