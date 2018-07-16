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

#if !defined(MPTOOLKIT_TENSOR_TENSOR_EIGEN_H)
#define MPTOOLKIT_TENSOR_TENSOR_EIGEN_H

#include "tensor.h"
#include "blas/diagonalmatrix.h"
#include "tensor_types.h"

// These are mostly unused - just some of them for infinite wavefunctions (maybe)

namespace Tensor
{

// singular value decomposition, singular values returned as a diagonal matrix.
// It is not necessary to initialize U,D,Vh with the correct dimensions,
// they will be move-assigned with new matrices anyway.
void
SVD(MatrixOperator const& m,
    MatrixOperator& U,
    RealDiagonalOperator& D,
    MatrixOperator& Vh);

// Version of the SVD where the dimension of D is the number of columns of m
void
SVD_FullColumns(MatrixOperator const& m,
		MatrixOperator& U,
		RealDiagonalOperator& D,
		MatrixOperator& Vh);

// Version of the SVD where the dimension of D is the number of rows of m
void
SVD_FullRows(MatrixOperator const& m,
	     MatrixOperator& U,
	     RealDiagonalOperator& D,
	     MatrixOperator& Vh);

// Version of the SVD where the dimension of D is max(m.rows(), m.cols())
void
SVD_Full(MatrixOperator const& m,
	 MatrixOperator& U,
	 RealDiagonalOperator& D,
	 MatrixOperator& Vh);

// SVD where every basis is a regular basis
void
SVD_Regular(MatrixOperator const& m,
	    MatrixOperator& U,
	    RealDiagonalOperator& D,
	    MatrixOperator& Vh);


// DiagonalizeHermitian
//
// U is modified to be the transform matrix, such that
// U' * Result' * herm(U') = U

RealDiagonalOperator
DiagonalizeHermitian(MatrixOperator& U);

// SqrtDiagonal
//
// Takes the square root of a real positive matrix.
// Elements that are < OrthoTol are assumed to be zero

RealDiagonalOperator
SqrtDiagonal(RealDiagonalOperator x, double OrthoTol);

// InvertDiagonal
//
// Takes the inverse root of a real positive matrix.
// Elements that are < OrthoTol are assumed to be zero

RealDiagonalOperator
InvertDiagonal(RealDiagonalOperator x, double OrthoTol);


#if 0

// DiagonalizeHermitian, diagonalizes 'x', returns the transform matrix
// such that x' = result' * x * herm(result') is diagonal.
// x must be a symmetric scalar operator, with a regular basis
//(ie. at most one subspace per quantum number).
// Note, this interface differs from LinearAlgebra::DiagonalizeSymmetric
template <typename T>
IrredTensor<T, VectorBasis, VectorBasis>
DiagonalizeHermitian(IrredTensor<T, VectorBasis, VectorBasis>& x);


LinearAlgebra::Vector<double>
EigenvaluesHermitian(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                     VectorBasis, VectorBasis> const& x);

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


// takes the square root of a positive diagonal matrix, with some allowance
// for a slightly negative entry (within -tol)
IrredTensor<LinearAlgebra::Matrix<std::complex<double> >, VectorBasis, VectorBasis>
SqrtDiagonal(IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
             VectorBasis, VectorBasis> const& Op, double Tol = 1E-15);

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

#endif

} // namespace Tensor

#include "tensor_eigen.cc"

#endif
