// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/regularize.h
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

#if !defined(REGULARIZE_H_DHCKJDHUREYT7845Y78Y78TY78TY78T)
#define REGULARIZE_H_DHCKJDHUREYT7845Y78Y78TY78TY78T

#include "tensor.h"
#include "blas/matrix.h"

namespace Tensor
{

// Regularizes the basis, returns the transform matrix.
// The regular basis is Result'.Basis1()
// Result'.Basis2() is b.
template <typename T>
IrredTensor<T, VectorBasis, VectorBasis>
Regularize(VectorBasis const& b);

// is_regular_basis
// Returns true if the basis is regular
bool is_regular_basis(VectorBasis const& b);

// regularize a BasisList
IrredTensor<blas::Matrix<double>, VectorBasis, BasisList>
Regularize(BasisList const& b);

// split a VectorBasis into a BasisList.  Returns the transpose of the tranform matrix
IrredTensor<blas::Matrix<double>, VectorBasis, BasisList>
SplitBasis(VectorBasis const& b);

// Map a SimpleOperator defined over a 1x1 matrix onto the element itself.
// This function is necessary because transforming a split basis into an operator
// over a BasisList doesn't result in an ordinary SimpleOperator, but something
// like a MatrixOperator but with 1x1 matrices.
// precondition: All elements of Op are 1x1 matrices.
template <typename T>
IrredTensor<T, BasisList, BasisList>
map_1x1_operator(IrredTensor<blas::Matrix<T>, BasisList, BasisList> const& Op);

   // convert a simple operator to a matrix operator
template <typename T>
IrredTensor<blas::Matrix<T>, VectorBasis, VectorBasis>
ToMatrixOperator(IrredTensor<T, BasisList, BasisList> const& Op);

} // namespace Tensor

#include "regularize.cc"

#endif
