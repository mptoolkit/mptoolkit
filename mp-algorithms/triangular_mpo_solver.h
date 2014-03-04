// -*- C++ -*- $Id$
//
// Solver for triangular MPO's

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
