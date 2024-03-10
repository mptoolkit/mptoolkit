// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// tensor/tensor_eigen.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// Return the singular values of a scalar tensor.  The singular values can be either the raw values
// (probably not useful)
// or they can be normalized by the quantum dimension of the quantum number sector, which is what you want
// if the normalization is with respect to the sum of the singular values,
// or they can be normalized by the square root of the quantum dimension, which is what you want if
// the normalization is with respect to the sum of the squares of the singular values.
enum class SingularValueNormalization { None, QDim, SqrtQDim };
LinearAlgebra::Vector<double>
SingularValues(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
               VectorBasis, VectorBasis> const& m, SingularValueNormalization n = SingularValueNormalization::None);

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

// QR factorization of a scalar operator, result is a pair of Q, R, where Q is rectangular, and R is upper-triangular.
// We require that the number of rows of A >= number of columns, in each sector. (FIXME: I think this isn't actually required)
// The Full version returns Q as m x m matrix, and R as m x n.
std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
QR_FactorizeFull(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> A);

// QR factorization of a scalar operator, result is a pair of Q, R, where Q is rectangular, and R is square upper-triangular.
// We require that the number of rows of A >= number of columns, in each sector.
// The default ('thin') version returns Q as an m x n matrix, and R as n x n.
// TODO: this ought to return Q as m x min(m,n)
std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
QR_Factorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> A);

// LQ factorization of a scalar operator, result is a pair of L, Q, where Q is rectangular, and L is lower-triangular.
// We require that the number of columns of A >= number of rows, in each sector.
// The Full version returns L as an m x n matrix, and Q as n x n.
std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
LQ_FactorizeFull(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> A);

// LQ factorization of a scalar operator, result is a pair of L, Q, where Q is rectangular, and L is lower-triangular.
// We require that the number of columns of A >= number of rows, in each sector.
// The default ('thin') version returns L as an m x m matrix, and Q as m x n.
// TODO: this ought to return Q as min(m,n) x n
std::pair<IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>,
         IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>>
LQ_Factorize(IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis> A);

// inverts a diagonal operator, with a given cutoff of singular values
IrredTensor<LinearAlgebra::Matrix<std::complex<double>>, VectorBasis, VectorBasis>
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
