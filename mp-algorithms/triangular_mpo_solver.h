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

#include "mps/momentum_operations.h"

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator.
// The LeftIdentity and RightIdentity are the right and left eigenmatrices of the 
// transfer operator in the Basis1 of Psi.
// If Psi is left-orthogonal then LeftIdentity = I and RightIdentity = Rho
// if Psi is right-orthogonal then LeftIdentity = rho and RightIdentity = I
KMatrixPolyType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularMPO const& Op, MatrixOperator const& LeftIdentity,
              MatrixOperator const& RightIdentity, int Verbose = 0);

// Solve an MPO in the right-handed sense, as Op * x_R = x_R * lambda
// We currently assume there is only one eigenvalue 1 of the transfer operator.
// The LeftIdentity and RightIdentity are the right and left eigenmatrices of the 
// transfer operator in the Basis2 of Psi.
// If Psi is left-orthogonal then LeftIdentity = I and RightIdentity = Rho
// if Psi is right-orthogonal then LeftIdentity = rho and RightIdentity = I
KMatrixPolyType
SolveMPO_Right(LinearWavefunction const& Psi, QuantumNumber const& QShift,
	       TriangularMPO const& Op, MatrixOperator const& LeftIdentity,
	       MatrixOperator const& RightIdentity, int Verbose = 0);

#endif
