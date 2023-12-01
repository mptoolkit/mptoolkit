// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/eigen.h
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
//
// Complex diagonalization wrapper by Thomas Barthel
//
// If a Range parameter is given, then only the eigenvalues/vectors from that range
// (ordered from smallest to largest) are calculated.
//


#if !defined(EIGEN_H_JDSHF849WYT879Y89RYT89YG89YU95UYGP)
#define EIGEN_H_JDSHF849WYT879Y89RYT89YG89YU95UYGP

#include "matrix.h"

namespace LinearAlgebra
{

//
// LinearSolveSPD
//
// Solves a symmetric, positive definite system of linear equations  M*v=rhs.
// For input M is the symmetric, positive def., real matrix
// and RHS is a matrix, with rows representing different right hand sides,
// for which the SoEq is to be solved.
// Returns the solutions v as column vectors.
// Postcondition:
// Satisfies m * result' = rhs

template <typename M1, typename M2,
          typename M1i = typename interface<M1>::type,
          typename M2i = typename interface<M2>::type>
struct ImplementLinearSolveSPD {};

template <typename M1, typename M2>
typename ImplementLinearSolveSPD<M1, M2>::result_type
LinearSolveSPD(M1 const& m, M2 const& rhs);

//
// LinearSolveHPD
//
// Solves a Hermitian, positive definite system of linear equations  M*v=rhs.
// For input M is the Hermitian, positive def., real matrix
// and RHS is a matrix, with rows representing different right hand sides,
// for which the SoEq is to be solved.
// On exit M contains the factor U or L from the Cholesky factorization
// and RHS the different solutions as rows.

template <typename M1, typename M2>
Matrix<std::complex<double>, ColMajor>
LinearSolveHPD(M1 const& m, M2 const& rhs,
               typename boost::enable_if<is_matrix<M1> >::type* = 0,
               typename boost::enable_if<is_matrix<M2> >::type* = 0);

//
// LinearSolve
//
// Solves a general system of linear equations  M*v=rhs.
// For input M is a generic matrix,
// and RHS is a matrix, with rows representing different right hand sides,
// for which the SoEq is to be solved.
// Returns the solutions v as column vectors.
// Postcondition:
// Satisfies m * result' = rhs

template <typename M1, typename M2>
Matrix<double, ColMajor>
LinearSolve(M1 const& m, M2 const& rhs,
            typename boost::enable_if<is_matrix<M1> >::type* = 0,
            typename boost::enable_if<is_matrix<M2> >::type* = 0);

//
// LeastSquares
//
// Least squares fit to a system of equations

template <typename M, typename V>
Vector<double>
LeastSquares(M const& m, V const& v,
            typename boost::enable_if<is_matrix<M> >::type* = 0,
            typename boost::enable_if<is_vector<V> >::type* = 0);

//
// LeastSquaresRegularized
//
// Least squares fit to a system of equations, with Tikhonov regularization.
// Returns the residual norm and the solution vector.
//
// Tikhonov regularization finds the solution x that minimizes ||Ax-b|| + alpha||x||
// where ||.|| is the vector 2-norm, and alpha is positive.
// This is implemented via SVD of A.
// Let A = UDV^\dagger, with D diagonal with entries d_i.
// Then the regularized solution is x = V X U^\dagger b,
// where X is diagonal with elements x_i = d_i / (d_i^2 + alpha).
//
// The residual is ||Ax-b|| + alpha||x|| = ||U DX U^\dagger x - b||
// The alpha||x|| term isn't included in the returned residual.
//

std::pair<double, Vector<std::complex<double>>>
LeastSquaresRegularized(Matrix<std::complex<double>> const& A, Vector<std::complex<double>> const& b,
                        double alpha = 1E-15);

//
// SingularValueDecomposition
//
// Singular values.  Returns the singular values in a LinearAlgebra::Vector<double>
//

template <typename A, typename D,
          typename Ai = typename interface<A>::type,
          typename Di = typename interface<D>::type>
struct ImplementSingularValues {};

template <typename A, typename D>
inline
typename ImplementSingularValues<A,D>::result_type
SingularValues(A const& a, D& d)
{
   return ImplementSingularValues<A,D>()(a,d);
}

// Singular value decomposition of a real matrix.  For input matrix A,
// D is a vector containing the singular values, and
// A = U * D * VT with U and VT near-orthogonal.

template <typename A, typename U, typename D, typename Vt,
          typename Ai = typename interface<A>::type,
          typename Ui = typename interface<U>::type,
          typename Di = typename interface<D>::type,
          typename Vti = typename interface<Vt>::type>
struct ImplementSingularValueDecomposition {};

template <typename A, typename U, typename D, typename Vt>
inline
typename ImplementSingularValueDecomposition<A,U,D,Vt>::result_type
SingularValueDecomposition(A const& a, U& u, D& d, Vt& vt)
{
   return ImplementSingularValueDecomposition<A,U,D,Vt>()(a,u,d,vt);
}

// Version that constructs only the left singular vectors

template <typename A, typename U, typename D,
          typename Ai = typename interface<A>::type,
          typename Ui = typename interface<U>::type,
          typename Di = typename interface<D>::type>
struct ImplementSingularValueDecompositionLeft {};

template <typename A, typename U, typename D>
inline
typename ImplementSingularValueDecompositionLeft<A,U,D>::result_type
SingularValueDecompositionLeft(A const& a, U& u, D& d)
{
   return ImplementSingularValueDecompositionLeft<A,U,D>()(a,u,d);
}

// Version that constructs only the right singular vectors

template <typename A, typename D, typename Vt,
          typename Ai = typename interface<A>::type,
          typename Di = typename interface<D>::type,
          typename Vti = typename interface<Vt>::type>
struct ImplementSingularValueDecompositionRight {};

template <typename A, typename D, typename Vt>
inline
typename ImplementSingularValueDecompositionRight<A,D,Vt>::result_type
SingularValueDecompositionRight(A const& a, D& d, Vt& vt)
{
   return ImplementSingularValueDecompositionRight<A,D,Vt>()(a,d,vt);
}

// Version that constructs the matrices to have dimension N x max(N,M) x M

template <typename A, typename U, typename D, typename Vt,
          typename Ai = typename interface<A>::type,
          typename Ui = typename interface<U>::type,
          typename Di = typename interface<D>::type,
          typename Vti = typename interface<Vt>::type>
struct ImplementSingularValueDecompositionFull {};

template <typename A, typename U, typename D, typename Vt>
inline
typename ImplementSingularValueDecompositionFull<A,U,D,Vt>::result_type
SingularValueDecompositionFull(A const& a, U& u, D& d, Vt& vt)
{
   return ImplementSingularValueDecompositionFull<A,U,D,Vt>()(a,u,d,vt);
}

// Given MxN matrix,
// U is MxM
// D is MxM
// Vt is MxN
template <typename A, typename U, typename D, typename Vt,
          typename Ai = typename interface<A>::type,
          typename Ui = typename interface<U>::type,
          typename Di = typename interface<D>::type,
          typename Vti = typename interface<Vt>::type>
struct ImplementSingularValueDecompositionFullLeft {};

template <typename A, typename U, typename D, typename Vt>
inline
typename ImplementSingularValueDecompositionFullLeft<A,U,D,Vt>::result_type
SingularValueDecompositionFullLeft(A const& a, U& u, D& d, Vt& vt)
{
   return ImplementSingularValueDecompositionFullLeft<A,U,D,Vt>()(a,u,d,vt);
}

// Given MxN matrix,
// U is MxN
// D is NxN
// Vt is NxN
template <typename A, typename U, typename D, typename Vt,
          typename Ai = typename interface<A>::type,
          typename Ui = typename interface<U>::type,
          typename Di = typename interface<D>::type,
          typename Vti = typename interface<Vt>::type>
struct ImplementSingularValueDecompositionFullRight {};

template <typename A, typename U, typename D, typename Vt>
inline
typename ImplementSingularValueDecompositionFullRight<A,U,D,Vt>::result_type
SingularValueDecompositionFullRight(A const& a, U& u, D& d, Vt& vt)
{
   return ImplementSingularValueDecompositionFullRight<A,U,D,Vt>()(a,u,d,vt);
}

// 'Full' singular value decomposition but only calculate the left vectors

template <typename A, typename U, typename D,
          typename Ai = typename interface<A>::type,
          typename Ui = typename interface<U>::type,
          typename Di = typename interface<D>::type>
struct ImplementSingularValueDecompositionLeftFull {};

template <typename A, typename U, typename D>
inline
typename ImplementSingularValueDecompositionLeftFull<A,U,D>::result_type
SingularValueDecompositionLeftFull(A const& a, U& u, D& d)
{
   return ImplementSingularValueDecompositionLeftFull<A,U,D>()(a,u,d);
}

// 'Full' singular value decomposition but only calculate the right vectors

template <typename A, typename D, typename Vt,
          typename Ai = typename interface<A>::type,
          typename Di = typename interface<D>::type,
          typename Vti = typename interface<Vt>::type>
struct ImplementSingularValueDecompositionRightFull {};

template <typename A, typename D, typename Vt>
inline
typename ImplementSingularValueDecompositionRightFull<A,D,Vt>::result_type
SingularValueDecompositionRightFull(A const& a, D& d, Vt& vt)
{
   return ImplementSingularValueDecompositionRightFull<A,D,Vt>()(a,d,vt);
}

//
// DiagonalizeSymmetric
//
// Diagonalizes a matrix in-place.  The matrix is replaced by the transform matrix, with
// the eigenvectors as sucessive row-vectors.  For input matrix M,
// X = M' is the transform matrix, and X * M * X^T is diagonal.
// Returns a vector of eigenvalues.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementDiagonalizeSymmetric {};

template <typename M>
typename ImplementDiagonalizeSymmetric<M>::result_type
DiagonalizeSymmetric(M& m);

template <typename M>
typename boost::enable_if<is_mutable_proxy<M>,
                          ImplementDiagonalizeSymmetric<M> >::type::result_type
DiagonalizeSymmetric(M const& m);

//
// DiagonalizeHermitian
//
// Diagonalizes a matrix in-place.  The matrix is replaced by the transform matrix, with
// the conjugate eigenvectors as sucessive row-vectors.  For input matrix M,
// X = M' is the transform matrix, and X * M * X^+ is diagonal.
// Returns a vector of eigenvalues.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementDiagonalizeHermitian {};

template <typename M>
typename ImplementDiagonalizeHermitian<M>::result_type
DiagonalizeHermitian(M& m);

template <typename M>
typename boost::enable_if<is_mutable_proxy<M>,
                          ImplementDiagonalizeHermitian<M> >::type::result_type
DiagonalizeHermitian(M const& m);

//
// Diagonalize
//
// Diagonalizes a matrix.  The left and right eigenvectors are returned
// as successive row vectors of Left and Right respectively.
//

template <typename M, typename L, typename R,
          typename Mi = typename interface<M>::type,
          typename Li = typename interface<L>::type,
          typename Ri = typename interface<R>::type>
struct ImplementDiagonalize {};

template <typename M, typename L, typename R>
typename ImplementDiagonalize<M, L&, R&>::result_type
Diagonalize(M const& m, L& Left, R& Right)
{
   return ImplementDiagonalize<M, L&, R&>::apply(m, Left, Right);
}

//
// EigenvaluesSymmetric, EigenvaluesHermitian
//
// returns the eigenvalues of a real-symmetric or complex-hermitian matrix
//

template <typename M>
Vector<double>
EigenvaluesSymmetric(M const& m);

template <typename M>
Vector<double>
EigenvaluesHermitian(M const& m);

//
// EigenvaluesComplex
//
// returns the eigenvalues of a general complex matrix
//

template <typename M>
Vector<std::complex<double>>
EigenvaluesComplex(M const& m);

//
// GeneralizedEigenSymmetric
//
// Solves the generalized eigenvalue problem A * x = lambda * B * x
// A and B  must be symmetric, and B must be positive definite.
// Matrices A and B are DESTROYED on exit (contain meaningless data)
//

template <typename A, typename B, typename Eigenval, typename Eigenvec,
          typename Ai = typename interface<A>::type,
          typename Bi = typename interface<B>::type,
          typename Eigenvali = typename interface<Eigenval>::type,
          typename Eigenveci = typename interface<Eigenvec>::type>
struct ImplementGeneralizedEigenSymmetric {};

template <typename A, typename B, typename Eigenval, typename Eigenvec>
inline
typename ImplementGeneralizedEigenSymmetric<A,B,Eigenval, Eigenvec>::result_type
GeneralizedEigenSymmetric(A const& a, B const& b,
                          Eigenval& eigenval, Eigenvec& eigenvec,
                          Range const& Which, double abstol = 2 * std::numeric_limits<double>::min())
{
   return ImplementGeneralizedEigenSymmetric<A,B,Eigenval, Eigenvec>()(a,b,eigenval, eigenvec, Which, abstol);
}

template <typename A, typename B, typename Eigenval, typename Eigenvec>
inline
typename ImplementGeneralizedEigenSymmetric<A,B,Eigenval, Eigenvec>::result_type
GeneralizedEigenHermitian(A const& a, B const& b,
                          Eigenval& eigenval, Eigenvec& eigenvec,
                          Range const& Which, double abstol = 2 * std::numeric_limits<double>::min())
{
   return ImplementGeneralizedEigenSymmetric<A,B,Eigenval, Eigenvec>()(a,b,eigenval, eigenvec, Which, abstol);
}

//
// CholeskyFactorize[Lower|Upper]
//
// Return the Cholesky factorization of a hermitian positive definite matrix.
// The factorization is in-place; the result vector is stored into the
// input matrix; only the upper or lower component is read and written to.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementCholeskyFactorizeUpper {};
template <typename M, typename Mi = typename interface<M>::type>
struct ImplementCholeskyFactorizeLower {};

template <typename M>
inline
typename ImplementCholeskyFactorizeUpper<M&>::result_type
CholeskyFactorizeUpper(M& m)
{
   return ImplementCholeskyFactorizeUpper<M&>()(m);
}

template <typename M>
inline
typename ImplementCholeskyFactorizeLower<M&>::result_type
CholeskyFactorizeLower(M& m)
{
   return ImplementCholeskyFactorizeLower<M&>()(m);
}

//
// SingularFactorize
//
// Return the factorization of a hermitian positive definite
// matrix into the 'square root' matrix X such that
// M = X * herm(X)
// via an eigenvalue decomposition.  This is more robust than
// the Cholesky factorization where the matrix is nearly
// singular.  Eigenvalues smaller than 0 are set to zero, so it is
// robust in the case of slightly non-positive matrices.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementSingularFactorize {};

template <typename M>
inline
typename ImplementSingularFactorize<M>::result_type
SingularFactorize(M const& m)
{
   return ImplementSingularFactorize<M>()(m);
}

//
// QR_FactorizeFull
//
// Performs the QR factorization of an MxN complex matrix with M >= N.
// The matrix is replaced by the upper triangular matrix R, and
// the unitary matrix Q is returned.  QR=X, where X is the original
// matrix, and R=X'.
// In the case where X is MxN and M > N, there is a choice for how to represent
// Q and R; we could choose full rank, Q is MxM, R is MxN, and is augmented by rows of zeros.
// Or we can choose minimal rank, Q is MxN, R is NxN.
// For the QR_FactorizeFull() function, we choose the first option, so that R is MxN and so doesn't need resizing.

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementQRFactorizeFull {};

template <typename M>
inline
typename ImplementQRFactorizeFull<M&>::result_type
QR_FactorizeFull(M& m)
{
   return ImplementQRFactorizeFull<M&>()(m);
}

// For the 'thin' QR factorization, we return a pair of matrices, and pass the original matrix by value, so that
// we can use move sematics to avoid creating a new matrix, where possible.  So we also abandon the generic interface
// and just write the functions directly.

std::tuple<Matrix<std::complex<double>>, Matrix<std::complex<double>>>
QR_Factorize(Matrix<std::complex<double>> M);

std::tuple<Matrix<double>, Matrix<double>>
QR_Factorize(Matrix<double> M);

//
// TridiagonalizeHermitian
//
// Tridiagonalizes a Hermitian matrix.  Input matrix is trashed on result.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementTridiagonalizeHermitian {};

template <typename M>
inline
typename ImplementTridiagonalizeHermitian<M&>::result_type
TridiagonalizeHermitian(M& m)
{
   return ImplementTridiagonalizeHermitian<M&>()(m);
}

//
// InvertHPD
//
// Inverts a Hermitian positive definite matrix in-place.
// Returns void.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementInvertHPD {};

template <typename M>
inline
typename ImplementInvertHPD<M&>::result_type
InvertHPD(M& m)
{
   return ImplementInvertHPD<M&>()(m);
}

//
// InvertGeneral
//
// Inverts a generic non-singular matrix in-place.
// Returns void.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementInvertGeneral {};

template <typename M>
inline
typename ImplementInvertGeneral<M&>::result_type
InvertGeneral(M& m)
{
   return ImplementInvertGeneral<M&>()(m);
}

//
// Invert[Lower|Upper]Triangular
//
// Inverts a upper or lower triangular matrix in-place.
// Returns void.
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementInvertLowerTriangular {};
template <typename M, typename Mi = typename interface<M>::type>
struct ImplementInvertUpperTriangular {};

template <typename M>
inline
typename ImplementInvertLowerTriangular<M&>::result_type
InvertLowerTriangular(M& m)
{
   return ImplementInvertLowerTriangular<M&>()(m);
}

template <typename M>
inline
typename ImplementInvertUpperTriangular<M&>::result_type
InvertUpperTriangular(M& m)
{
   return ImplementInvertUpperTriangular<M&>()(m);
}

#if defined(CONFIG_EXPOKIT)

#endif //  defined(CONFIG_EXPOKIT)

} // namespace LinearAlgebra

#include "eigen.cc"

namespace LinearAlgebra
{

template <typename M>
typename ImplementDiagonalizeSymmetric<M>::result_type
DiagonalizeSymmetric(M& m)
{
   return ImplementDiagonalizeSymmetric<M>::apply(m);
}

template <typename M>
typename boost::enable_if<is_mutable_proxy<M>,
                          ImplementDiagonalizeSymmetric<M> >::type::result_type
DiagonalizeSymmetric(M const& m)
{
   return ImplementDiagonalizeSymmetric<M>::apply(const_cast<M&>(m));
}

template <typename M1, typename M2>
typename ImplementLinearSolveSPD<M1, M2>::result_type
LinearSolveSPD(M1 const& m, M2 const& rhs)
{
   return ImplementLinearSolveSPD<M1, M2>::apply(m, rhs);
}

template <typename M>
typename ImplementDiagonalizeHermitian<M>::result_type
DiagonalizeHermitian(M& m)
{
   return ImplementDiagonalizeHermitian<M>::apply(m);
}

template <typename M>
typename boost::enable_if<is_mutable_proxy<M>,
                          ImplementDiagonalizeHermitian<M> >::type::result_type
DiagonalizeHermitian(M const& m)
{
   return ImplementDiagonalizeHermitian<M>::apply(const_cast<M&>(m));
}

} // namespace LinearAlgebra

#endif
