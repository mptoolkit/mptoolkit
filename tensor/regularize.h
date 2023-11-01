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

namespace Tensor
{

// Class to keep track of the index mapping for mapping a VectorBasis onto a regular basis
class Regularizer
{
   public:
      Regularizer() = delete;

      explicit Regularizer(VectorBasis const& b);

      // returns the index into the regular basis of the component i
      int IndexOf(int i) const;

      // returns the range mapping of the component i in the regular basis
      LinearAlgebra::Range RangeOf(int i) const;

      VectorBasis const& Basis() const { return RegularBasis; }

      VectorBasis const& OriginalBasis() const { return IrregularBasis; }


   private:
      std::vector<int>                  BasisMappingIndex;
      std::vector<LinearAlgebra::Range> BasisMappingRange;
      VectorBasis                       IrregularBasis;
      VectorBasis                       RegularBasis;
};

// Regularize Basis1 of a tensor, more efficiently than multiplying by the transformation tensor.
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis1(Regularizer const& R, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M);

// Version where we construct the regularizer implicitly
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis1(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M);

// Regularize Basis2 of a tensor, more efficiently than multiplying by the transformation tensor.
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis2(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R);

// Version where we construct the regularizer implicitly
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis2(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M);

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
RegularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R2);

// Not yet implemented
template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
RegularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& M, Regularizer const& R2);

// Inverse of Regularize
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
UnregularizeBasis1(Regularizer const& R, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M);

// Regularize Basis2 of a tensor, more efficiently than multiplying by the transformation tensor.
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
UnregularizeBasis2(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R);

template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
UnregularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> const& M, Regularizer const& R2);

// Not yet implemented
template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
UnregularizeBasis12(Regularizer const& R1, IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& M, Regularizer const& R2);

// Legacy Regularize for BasisList.  This would be better viewed as mapping a BasisList to a VectorBasis,
// with the inverse operation being SplitBasis (better renamed as something else?)

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList>
Regularize(BasisList const& b);

// is_regular_basis
// Returns true if the basis is regular
bool is_regular_basis(VectorBasis const& b);

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

   // convert a simple operator to a matrix operator
template <typename T>
IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis>
ToMatrixOperator(IrredTensor<T, BasisList, BasisList> const& Op);

} // namespace Tensor

#include "regularize.cc"

#endif
