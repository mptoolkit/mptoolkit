// -*- C++ -*- $Id$
//
// Solver for triangular MPO's

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TRIANGULAR_MPO_SOLVER_H)
#define MPTOOLKIT_MP_ALGORITHMS_TRIANGULAR_MPO_SOLVER_H

#include "mps/momentum_operations.h"

// Solve an MPO in the left-handed sense, as x_L * Op = lambda * x_L
// We currently assume there is only one eigenvalue 1 of the transfer operator
KMatrixPolyType
SolveMPO_Left(LinearWavefunction const& Psi, QuantumNumber const& QShift,
              TriangularMPO const& Op, MatrixOperator const& Rho,
              MatrixOperator const& Identity, int Verbose = 0);


#endif
