// -*- C++ -*- $Id$

#if !defined(REGULARIZE_H_DHCKJDHUREYT7845Y78Y78TY78TY78T)
#define REGULARIZE_H_DHCKJDHUREYT7845Y78Y78TY78TY78T

#include "tensor.h"

namespace Tensor
{

// Regularizes the basis, returns the transform matrix.
// The regular basis is Result'.Basis1()
// Result'.Basis2() is b.
IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
Regularize(VectorBasis const& b);

// is_regular_basis
// Returns true if the basis is regular
bool is_regular_basis(VectorBasis const& b);

// regularize a BasisList
IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList>
Regularize(BasisList const& b);

// split a VectorBasis into a BasisList.  Returns the transpose of the tranform matrix
IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList>
SplitBasis(VectorBasis const& b);

// Map a SimpleOperator defined over a 1x1 matrix onto the element itself.
// This function is necessary because transforming a split basis into an operator
// over a BasisList doesn't result in an ordinary SimpleOperator, but something
// like a MatrixOperator but with 1x1 matrices.
// precondition: All elements of Op are 1x1 matrices.
template <typename T>
IrredTensor<T, BasisList, BasisList>
map_1x1_operator(IrredTensor<LinearAlgebra::Matrix<T>, BasisList, BasisList> const& Op);

} // namespace Tensor

#include "regularize.cc"

#endif
