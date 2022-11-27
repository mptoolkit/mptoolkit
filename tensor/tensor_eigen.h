// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensor_eigen.h
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

#if !defined(TENSOR_EIGEN_H_DHCKJDHUREYT7845Y78Y78TY78TY78T)
#define TENSOR_EIGEN_H_DHCKJDHUREYT7845Y78Y78TY78TY78T

#include "tensor.h"
#include "linearalgebra/diagonalmatrix.h"

namespace Tensor
{

// DiagonalizeHermitian, diagonalizes 'x', returns the transform matrix
// such that x' = result' * x * herm(result') is diagonal.
// x must be a symmetric scalar operator, with a regular basis
//(ie. at most one subspace per quantum number).
// Note, this interface differs from LinearAlgebra::DiagonalizeSymmetric
IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
DiagonalizeHermitianInPlace(IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>& x);

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
DiagonalizeHermitianInPlace(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                     VectorBasis, VectorBasis>& x);

// This version returns a diagonal matrix
std::tuple<IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>,
           IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>>
DiagonalizeHermitian(IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis> x);

std::tuple<IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>,
           IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
DiagonalizeHermitian(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> x);

void
InvertHPD(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x);

void
InvertGeneral(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>& x);

void InvertIrregularHPD(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                        VectorBasis, VectorBasis>& x);

LinearAlgebra::Vector<double>
SingularValues(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& m);

void
SingularValueDecomposition(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& Vh);

// version to use explicitly a diagonal matrix
void
SingularValueDecomposition(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis> const& m,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& U,
                           IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                           IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                       VectorBasis, VectorBasis>& Vh);

// Version that keeps all of the zero singular values
void
SingularValueDecompositionFull(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                   VectorBasis, VectorBasis> const& m,
                               IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                           VectorBasis, VectorBasis>& U,
                               IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                               IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                           VectorBasis, VectorBasis>& Vh);

// Keep the Basis1 intact, that is, given input MxN matrix,
// U is MxM
// D is MxM
// Vh is MxN
void
SingularValueDecompositionKeepBasis1(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                   VectorBasis, VectorBasis> const& m,
                               IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                           VectorBasis, VectorBasis>& U,
                               IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                               IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                           VectorBasis, VectorBasis>& Vh);

// Keep the Basis2 intact, that is, given input MxN matrix,
// U is MxN
// D is NxN
// Vh is NxN
void
SingularValueDecompositionKeepBasis2(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                   VectorBasis, VectorBasis> const& m,
                               IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                           VectorBasis, VectorBasis>& U,
                               IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
                                           VectorBasis, VectorBasis, DiagonalStructure>& D,
                               IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                           VectorBasis, VectorBasis>& Vh);

// If the basis is already regular, we can avoid some book-keeping
void
SingularValueDecompositionRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                  VectorBasis, VectorBasis> const& m,
                                  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                  VectorBasis, VectorBasis>& U,
                                     IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                                 VectorBasis, VectorBasis>& D,
                                  IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                                              VectorBasis, VectorBasis>& Vh);

// SingularValueDecomposition of a SimpleOperator
void
SingularValueDecomposition(IrredTensor<std::complex<double>, BasisList, BasisList> const& m,
                           IrredTensor<std::complex<double>, BasisList, BasisList>& U,
                           IrredTensor<std::complex<double>, BasisList, BasisList>& D,
                           IrredTensor<std::complex<double>, BasisList, BasisList>& Vh);

// inverts a diagonal operator, with a given cutoff of singular values
IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& Op, double Cutoff);

// inverts a diagonal operator, with a default cutoff of singular values
IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& Op);

// Inversion of an explicitly diagonal matrix
template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<T>,
               VectorBasis, VectorBasis> const& Op, double Cutoff);

template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis>
InvertDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<T>,
               VectorBasis, VectorBasis> const& Op);

// Take the diagonal part of a tensor
template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
ExtractDiagonal(IrredTensor<LinearAlgebra::Matrix<T>, VectorBasis, VectorBasis> A);

// Take the real diagonal part of a tensor
template <typename T>
IrredTensor<LinearAlgebra::DiagonalMatrix<T>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
ExtractDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<T>>, VectorBasis, VectorBasis> A);

// takes the square root of a positive diagonal matrix, with some allowance
// for a slightly negative entry (within -tol)
IrredTensor<LinearAlgebra::DiagonalMatrix<std::complex<double>>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
SqrtDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<std::complex<double>>,
             VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& Op, double Tol = 1E-14);

IrredTensor<LinearAlgebra::DiagonalMatrix<double>, VectorBasis, VectorBasis, Tensor::DiagonalStructure>
SqrtDiagonal(IrredTensor<LinearAlgebra::DiagonalMatrix<double>,
             VectorBasis, VectorBasis, Tensor::DiagonalStructure> const& Op, double Tol = 1E-14);

// Cholesky factorization of a hermitian positive definite matrix.  We assume here that
// Basis1 == Basis2

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeUpperRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                              VectorBasis, VectorBasis> const& m);

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
CholeskyFactorizeUpper(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                       VectorBasis, VectorBasis> const& m);

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SingularFactorizeRegular(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                         VectorBasis, VectorBasis> const& m);

IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SingularFactorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                  VectorBasis, VectorBasis> const& m);

} // namespace Tensor

#include "tensor_eigen.cc"

#endif
