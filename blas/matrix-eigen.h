// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-lapack.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_MATRIX_EIGEN_H)
#define MPTOOLKIT_BLAS_MATRIX_EIGEN_H

#include <complex>
#include "detail/matrix-lapack.h"

namespace blas
{

//
// middle-layer LAPACK wrappers, that forward from a matrix/vector ref to low-level storage
//

//
// Diagonalize
//
// Diagonalizes a matrix.  The left and right eigenvectors are returned
// as successive column vectors of Left and Right respectively.
// The input matrix is destroyed, so we have two versions were we either copy
// the input matrix, or if it is an rvalue-reference, we can destroy it.
/// On input, M is an NxN matrix, V is an N-dimensional vector,
// and Left and Right are NxN matrices.
//
// TODO: add NormalVectorProxy versions

template <typename U, typename V, typename W, typename X, typename Tag>
inline
void Diagonalize(NormalMatrix<std::complex<double>, U, Tag>&& M,
		 NormalVector<std::complex<double>, V, Tag>& v,
		 NormalMatrix<std::complex<double>, W, Tag>& Left,
		 NormalMatrix<std::complex<double>, X, Tag>& Right)
{
   CHECK_EQUAL(M.rows(), M.cols());
   CHECK_EQUAL(M.rows(), v.size());
   CHECK_EQUAL(M.rows(), Left.rows());
   CHECK_EQUAL(M.rows(), Left.cols());
   CHECK_EQUAL(M.rows(), Right.rows());
   CHECK_EQUAL(M.rows(), Right.cols());
   detail::Diagonalize(M.rows(), std::move(M).storage(), M.leading_dimension(), v.storage(),
		       Left.storage(), Left.leading_dimension(),
		       Right.storage(), Right.leading_dimension());
}

template <typename U, typename V, typename W, typename X, typename Tag>
inline
void Diagonalize(NormalMatrix<std::complex<double>, U, Tag> const& M,
		 NormalVector<std::complex<double>, V, Tag>& v,
		 NormalMatrix<std::complex<double>, W, Tag>& Left,
		 NormalMatrix<std::complex<double>, X, Tag>& Right)
{
   Diagonalize(copy(M), v, Left, Right);
}

//
// DiagonalizeSymmetric
//
// Diagonalizes a matrix in-place.  The matrix is replaced by the transform matrix, with
// the eigenvectors as sucessive column-vectors.  For input matrix M,
// X = M' is the transform matrix, E is the diagonal matrix of eigenvalues,
// we have MX = XE
//
// The vector v must be initialized to the correct size.

template <typename U, typename V, typename Tag>
inline
void DiagonalizeSymmetric(NormalMatrix<double, U, Tag>& M, NormalVector<double, V, Tag>& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   detail::DiagonalizeSymmetric(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

// Version that takes a proxy-reference for the eigenvalues
template <typename U, typename V, typename Tag>
inline
void DiagonalizeSymmetric(NormalMatrix<double, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   detail::DiagonalizeSymmetric(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

// TODO: we could also add versions where M is passed as a NormalMatrixProxy

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<double, U, Tag>& M, NormalVector<double, V, Tag>& v)
{
   DiagonalizeSymmetric(M.as_derived(), v.as_derived());
}

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<double, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   DiagonalizeSymmetric(M.as_derived(), static_cast<V&&>(v.as_derived()));
}

//
// complex
//

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<std::complex<double>, U, Tag>& M, NormalVector<double, V, Tag>& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   TRACE(M.as_derived())(M.rows())(M.leading_dimension());
   detail::DiagonalizeHermitian(M.rows(), M.storage(), M.leading_dimension(), v.storage());
}

template <typename U, typename V, typename Tag>
inline
void DiagonalizeHermitian(NormalMatrix<std::complex<double>, U, Tag>& M, NormalVectorProxy<double, V, Tag>&& v)
{
   CHECK_EQUAL(M.rows(), M.cols());
   detail::DiagonalizeHermitian(M.rows(), M.storage(), M.leading_dimension(), std::move(v).storage());
}

//
// singular value decomposition
//
// The input matrix is destroyed

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V, typename Tag>
void
SVD(NormalMatrix<Scalar, M, Tag>&& Mmat, NormalMatrix<Scalar, U, Tag>& Umat,
    NormalVector<Real, D, Tag>& Dvec, NormalMatrix<Scalar, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), std::min(Mmat.rows(), Mmat.cols()));
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
				      Umat.storage(), Umat.leading_dimension(),
				      Vmat.storage(), Vmat.leading_dimension());
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V, typename Tag>
void
SVD(NormalMatrix<Scalar, M, Tag>&& Mmat, NormalMatrix<Scalar, U, Tag>& Umat,
    NormalVectorProxy<Real, D, Tag>&& Dvec, NormalMatrix<Scalar, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), std::min(Mmat.rows(), Mmat.cols()));
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.rows(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
				      std::move(Dvec).storage(),
				      Umat.storage(), Umat.leading_dimension(),
				      Vmat.storage(), Vmat.leading_dimension());
}


// SVD_FullRows
//
// Given M as an m*n matrix, U is m*m, D is m*m, V is m*n
//

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V, typename Tag>
void
SVD_FullRows(NormalMatrix<Scalar, M, Tag>&& Mmat, NormalMatrix<Scalar, U, Tag>& Umat,
	     NormalVector<Real, D, Tag>& Dvec, NormalMatrix<Scalar, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.rows());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.rows() < Mmat.cols())
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
					 Umat.storage(), Umat.leading_dimension(),
					 Vmat.storage(), Vmat.leading_dimension());
   }
   else
   {
      // need to zero the additional elements in Dvec
      clear(Dvec);
      detail::SingularValueDecompositionFull(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
					     Umat.storage(), Umat.leading_dimension(),
					     Vmat.storage(), Vmat.leading_dimension());
   }
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V, typename Tag>
void
SVD_FullRows(NormalMatrix<Scalar, M, Tag>&& Mmat, NormalMatrix<Scalar, U, Tag>& Umat,
	     NormalVectorProxy<Real, D, Tag>&& Dvec, NormalMatrix<Scalar, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.rows());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.rows() < Mmat.cols())
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
					 std::move(Dvec).storage(),
					 Umat.storage(), Umat.leading_dimension(),
					 Vmat.storage(), Vmat.leading_dimension());
   }
   else
   {
      // need to zero the additional elements in Dvec
      clear(Dvec);
      detail::SingularValueDecompositionFull(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
					     std::move(Dvec).storage(),
					     Umat.storage(), Umat.leading_dimension(),
					     Vmat.storage(), Vmat.leading_dimension());
   }
}

// SVD_FullCols
//
// Given M as an m*n matrix, U is m*n, D is n*n, V is n*n
//

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V, typename Tag>
void
SVD_FullCols(NormalMatrix<Scalar, M, Tag>&& Mmat, NormalMatrix<Scalar, U, Tag>& Umat,
	     NormalVector<Real, D, Tag>& Dvec, NormalMatrix<Scalar, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.cols());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   if (Mmat.cols() < Mmat.rows())
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(),
                                         Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
					 Umat.storage(), Umat.leading_dimension(),
					 Vmat.storage(), Vmat.leading_dimension());
   }
   else
   {
      // need to zero the additional elements in Dvec
      clear(Dvec);
      detail::SingularValueDecompositionFull(Mmat.rows(), Mmat.cols(),
                                             Mmat.storage(), Mmat.leading_dimension(), Dvec.storage(),
					     Umat.storage(), Umat.leading_dimension(),
					     Vmat.storage(), Vmat.leading_dimension());
   }
}

template <typename Scalar, typename Real, typename M, typename U, typename D, typename V, typename Tag>
void
SVD_FullCols(NormalMatrix<Scalar, M, Tag>&& Mmat, NormalMatrix<Scalar, U, Tag>& Umat,
	     NormalVectorProxy<Real, D, Tag>&& Dvec, NormalMatrix<Scalar, V, Tag>& Vmat)
{
   CHECK_EQUAL(Dvec.size(), Mmat.cols());
   CHECK_EQUAL(Mmat.rows(), Umat.rows());
   CHECK_EQUAL(Umat.cols(), Dvec.size());
   CHECK_EQUAL(Dvec.size(), Vmat.rows());
   CHECK_EQUAL(Vmat.cols(), Mmat.cols());
   //TRACE("WWWW")(Mmat);
   if (Mmat.cols() < Mmat.rows())
   {
      detail::SingularValueDecomposition(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
					 std::move(Dvec).storage(),
					 Umat.storage(), Umat.leading_dimension(),
					 Vmat.storage(), Vmat.leading_dimension());
   }
   else
   {
      // need to zero the additional elements in Dvec
      clear(Dvec);
      detail::SingularValueDecompositionFull(Mmat.rows(), Mmat.cols(), Mmat.storage(), Mmat.leading_dimension(),
					     std::move(Dvec).storage(),
					     Umat.storage(), Umat.leading_dimension(),
					     Vmat.storage(), Vmat.leading_dimension());
   }
   TRACE(Umat)(Vmat);
}

//
// LinearSolve
//
// Solves a general system of linear equations  M*v=rhs.
// For input M is a general square matrix,
// and RHS is a matrix, with rows representing different right hand sides
// for which the system of equations is to be solved.
// On exit, the RHS matrix is replaced by the solution vectors.
// Satisfies m * result' = rhs

// version with one RHS vector
template <typename U, typename V, typename Tag>
void
LinearSolve(NormalMatrix<double, U, Tag> const& M, NormalVector<double, V, Tag>& RHS);

// version with multiple RHS's
template <typename U, typename V, typename Tag>
void
LinearSolve(NormalMatrix<double, U, Tag> const& M, NormalMatrix<double, V, Tag>& RHS);

template <typename U, typename V, typename Tag>
void
LinearSolve(NormalMatrix<std::complex<double>, U, Tag> const& M,
	    NormalVector<std::complex<double>, V, Tag>& RHS);

// version with multiple RHS's
template <typename U, typename V, typename Tag>
void
LinearSolve(NormalMatrix<std::complex<double>, U, Tag> const& M,
	    NormalMatrix<std::complex<double>, V, Tag>& RHS);


} // namespace blas


#endif
