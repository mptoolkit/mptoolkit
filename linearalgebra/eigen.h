// -*- C++ -*- $Id$
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
LinearSolveSPD(M1 const& m, M2 const& rhs)
{
   return ImplementLinearSolveSPD<M1, M2>::apply(m, rhs);
}

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
// SingularValueDecomposition
//
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
// QR_Factorize
//
// Performs the QR factorization of a complex matrix.
// The matrix is replaced by the upper triangular matrix R, and
// the unitary matrix Q is returned.  QR=X, where X is the original
// matrix, and R=X'
//

template <typename M, typename Mi = typename interface<M>::type>
struct ImplementQRFactorize {};

template <typename M>
inline
typename ImplementQRFactorize<M&>::result_type
QR_Factorize(M& m)
{
   return ImplementQRFactorize<M&>()(m);
}

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

#endif
