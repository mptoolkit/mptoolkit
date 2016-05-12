// -*- C++ -*- $Id$
//
// Solver for triangular MPO's
//
// Assumptions:
// We assume that the identity is in the scalar sector of the transfer matrix.
// This should be true because we form the transfer matrix of a single wavefunction,
// we don't allow <psi1|X|psi2>.
// We assume that there is at most one eigenvalue 1 in each symmetry sector.  This allows
// for wavefunctions that are non-injective due to symmetry (eg Z2 symmetry in the ordered phase).
// We also allow for MPO's that have a non-trivial unitary on the diagonal, but again it should
// have at most one eigenvalue 1 per symmetry sector.
//


#if !defined(MPTOOLKIT_MP_ALGORITHMS_TRIANGULAR_MPO_SOLVER_H)
#define MPTOOLKIT_MP_ALGORITHMS_TRIANGULAR_MPO_SOLVER_H

#include "wavefunction/momentum_operations.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/infinitewavefunctionright.h"

// For identifying an eigenvalue 1 of the transfer matrix, we need
// an epsilon tolerance.  1E-12 proved to be a bit too small in some
// cases.
double const DefaultEigenUnityEpsilon = 1E-12;

// default tolerance for eigensolver and linear solver.  This is also a bit small in some cases.
double const DefaultTol = 1E-14;

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator.
// The LeftIdentity and RightIdentity are the right and left eigenmatrices of the 
// transfer operator in the Basis1 of Psi.
// If Psi is left-orthogonal then LeftIdentity = I and RightIdentity = Rho
// if Psi is right-orthogonal then LeftIdentity = rho and RightIdentity = I
// The return value is the complete set of E matrices.  The final expectation value
// is contained in the last element.
// If NeedFinalMatrix is false, then we don't calculate the final matrix elements
// (components perpendicular to the identity) for the last component.  This corresponds
// to the case where we want to calculate an expectation value only, and we
// don't need the complete matrix elements.
// The EMatK is a vector of the already-converged E-matrix elements.
// If no matrix elements are known then this can be initialized to the empty matrix.
// If EMatK is non-empty, then the supplied elements MUST be exactly the 
// elements of the final solution (this isn't checked, but the behaviour 
// will be unspecified otherwise).
// This is intended, eg for calculating higher powers of an MPO, where the solution at a lower
// power can be re-used to speed up the next power.
// The final matrix element is EMatK.back()' on exit, the expectation value is the
// overlap of this matrix element with the density matrix (RightIdentity).
void
SolveMPO_Left(std::vector<KMatrixPolyType>& EMatK,
	      LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularMPO const& Op, MatrixOperator const& LeftIdentity,
              MatrixOperator const& RightIdentity, 
	      bool NeedFinalMatrix, double Tol = DefaultTol,
	      double EigenUnityEpsilon = DefaultEigenUnityEpsilon, int Verbose = 0);

//
// 'Simple' solvers, for first order operators (eg Hamiltonians).  Return value is
// the eigenvalue per wavefunction unit cell.  The wavefunction must be in the appropriate
// left or right canonical form.  Prefer the versions that take an InfiniteWavefunctionLeft/Right
//

std::complex<double>
SolveSimpleMPO_Left(StateComponent& E, LinearWavefunction const& Psi,
		    QuantumNumber const& QShift, TriangularMPO const& Op,
		    MatrixOperator const& Rho, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveSimpleMPO_Left(StateComponent& E, InfiniteWavefunctionLeft const& Psi,
		    TriangularMPO const& Op, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveSimpleMPO_Right(StateComponent& F, LinearWavefunction const& Psi,
		     QuantumNumber const& QShift, TriangularMPO const& Op,
		     MatrixOperator const& Rho, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveSimpleMPO_Right(StateComponent& F, InfiniteWavefunctionRight const& Psi,
		     TriangularMPO const& Op, double Tol = DefaultTol, int Verbose = 0);

#endif
