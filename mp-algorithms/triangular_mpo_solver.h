// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/triangular_mpo_solver.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// If Degree is set, then set the maximum degree to this, and cull elements of a larger
// degree (with a warning if the matrix element is numerically non-zero).
void
SolveMPO_Left(std::vector<KMatrixPolyType>& EMatK,
              LinearWavefunction const& Psi, QuantumNumber const& QShift,
              BasicTriangularMPO const& Op, MatrixOperator const& Identity,
              MatrixOperator const& Rho,
              bool NeedFinalMatrix, int Degree = 0, double Tol = DefaultTol,
              double EigenUnityEpsilon = DefaultEigenUnityEpsilon, int Verbose = 0);

// Solve the cross expectation value <psi1|Op|psi2> / <psi1|psi2>.
// LeftIdentity and RightIdentity are a left/right eigenpair of the mixed transfer matrix,
// in the Basis1() basis of psi1/psi2.  lambda is the transfer matrix eigenvalue.
void
SolveMPO_Left_Cross(std::vector<KMatrixPolyType>& EMatK,
                    LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                    BasicTriangularMPO const& Op, MatrixOperator const& Identity,
                    MatrixOperator const& Rho,
                    std::complex<double> lambda,
                    bool NeedFinalMatrix, int Degree = 0, double Tol = DefaultTol,
                    double EigenUnityEpsilon = DefaultEigenUnityEpsilon, int Verbose = 0);

// MPO solver for 'simple' triangular MPO's.  These are a bit mis-named,
// they are not so simple, but refers to an MPO that has no momentum dependence,
// and if there are any diagonal components they are proportional to the identity,
// with either prefactor 1.0 or |prefactor| < 1.0

void
SolveSimpleMPO_Left(std::vector<MatrixPolyType>& EMat,
                   LinearWavefunction const& Psi, QuantumNumber const& QShift,
                   BasicTriangularMPO const& Op,
                   MatrixOperator const& Identity,
                   MatrixOperator const& Rho, bool NeedFinalMatrix,
                   int Degree, double Tol,
                   double UnityEpsilon, int Verbose);

void
SolveSimpleMPO_Right(std::vector<MatrixPolyType>& FMat,
                   LinearWavefunction const& Psi, QuantumNumber const& QShift,
                   BasicTriangularMPO const& Op,
                   MatrixOperator const& Identity,
                   MatrixOperator const& Rho, bool NeedFinalMatrix,
                   int Degree, double Tol,
                   double UnityEpsilon, int Verbose);

// Hamiltonian operators.  These are not necessarily strictly first order, so
// it is implemented as the fixed point polynomial evaluated at n=0.

std::complex<double>
SolveHamiltonianMPO_Left(StateComponent& E, LinearWavefunction const& Psi,
                       QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                       MatrixOperator const& Rho, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveHamiltonianMPO_Left(StateComponent& E, InfiniteWavefunctionLeft const& Psi,
                   BasicTriangularMPO const& Op, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveHamiltonianMPO_Right(StateComponent& F, LinearWavefunction const& Psi,
                    QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                    MatrixOperator const& Rho, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveHamiltonianMPO_Right(StateComponent& F, InfiniteWavefunctionRight const& Psi,
                    BasicTriangularMPO const& Op, double Tol = DefaultTol, int Verbose = 0);

//
// Example code for first order operators.  These work for almost all Hamiltonian operators,
// although in practice we use the SimpleMPO variants above since some exotic cases (eg gauge fields)
// end up having a higher order MPO, although the final expectation value is linear.
// This code is left for demonstration purposes.
//

std::complex<double>
SolveFirstOrderMPO_Left(StateComponent& E, LinearWavefunction const& Psi,
                        QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                        MatrixOperator const& Rho, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveFirstOrderMPO_Left(StateComponent& E, InfiniteWavefunctionLeft const& Psi,
                    BasicTriangularMPO const& Op, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveFirstOrderMPO_Right(StateComponent& F, LinearWavefunction const& Psi,
                     QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                     MatrixOperator const& Rho, double Tol = DefaultTol, int Verbose = 0);

std::complex<double>
SolveFirstOrderMPO_Right(StateComponent& F, InfiniteWavefunctionRight const& Psi,
                     BasicTriangularMPO const& Op, double Tol = DefaultTol, int Verbose = 0);

// Solvers for excitation ansatz wavefunctions, where PsiTri is the
// "triangular" unit cell and ExpIK is the complex phase per MPS unit cell.

void
SolveMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK, std::vector<KMatrixPolyType> const& CTriK,
                 LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                 QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                 MatrixOperator const& TLeft, MatrixOperator const& TRight,
                 std::complex<double> ExpIK, int Degree, double Tol = DefaultTol, double UnityEpsilon = DefaultEigenUnityEpsilon,
                 bool NeedFinalMatrix = true, bool EAOptimization = false, int Verbose = 0);

void
SolveMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK, std::vector<KMatrixPolyType> const& CTriK,
                  LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                  QuantumNumber const& QShift, BasicTriangularMPO const& Op,
                  MatrixOperator const& TLeft, MatrixOperator const& TRight,
                  std::complex<double> ExpIK, int Degree, double Tol = DefaultTol, double UnityEpsilon = DefaultEigenUnityEpsilon,
                  bool NeedFinalMatrix = true, bool EAOptimization = false, int Verbose = 0);

void
SolveMPO_EA_Left(std::vector<KMatrixPolyType>& EMatK, std::vector<KMatrixPolyType> const& CTriK,
                 LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                 QuantumNumber const& QShift, ProductMPO const& Op,
                 MatrixOperator const& TLeft, MatrixOperator const& TRight,
                 std::complex<double> ExpIK, int Degree, double Tol = DefaultTol, double UnityEpsilon = DefaultEigenUnityEpsilon,
                 bool NeedFinalMatrix = true, bool EAOptimization = false, int Verbose = 0);

void
SolveMPO_EA_Right(std::vector<KMatrixPolyType>& FMatK, std::vector<KMatrixPolyType> const& CTriK,
                  LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                  QuantumNumber const& QShift, ProductMPO const& Op,
                  MatrixOperator const& TLeft, MatrixOperator const& TRight,
                  std::complex<double> ExpIK, int Degree, double Tol = DefaultTol, double UnityEpsilon = DefaultEigenUnityEpsilon,
                  bool NeedFinalMatrix = true, bool EAOptimization = false, int Verbose = 0);

std::vector<KMatrixPolyType>
CalculateCTriK_Left(std::vector<KMatrixPolyType> const& EMatKNorth, std::vector<KMatrixPolyType> const& EMatKEast,
                    std::vector<KMatrixPolyType> const& EMatKNorthEast, LinearWavefunction const& Psi1,
                    LinearWavefunction const& Psi2, LinearWavefunction const& PsiTri1,
                    LinearWavefunction const& PsiTri2, QuantumNumber const& QShift,
                    GenericMPO const& Op, std::complex<double> ExpIK1, std::complex<double> ExpIK2);

std::vector<KMatrixPolyType>
CalculateCTriK_Right(std::vector<KMatrixPolyType> const& FMatKSouth, std::vector<KMatrixPolyType> const& FMatKWest,
                     std::vector<KMatrixPolyType> const& FMatKSouthWest, LinearWavefunction const& Psi1,
                     LinearWavefunction const& Psi2, LinearWavefunction const& PsiTri1,
                     LinearWavefunction const& PsiTri2, QuantumNumber const& QShift,
                     GenericMPO const& Op, std::complex<double> ExpIK1, std::complex<double> ExpIK2);

void
SolveHamiltonianMPO_EA_Left(StateComponent& E1, StateComponent const& E0,
                            LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                            LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                            BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                            MatrixOperator const& TRight, std::complex<double> ExpIK,
                            double Tol = DefaultTol, int Verbose = 0);

void
SolveHamiltonianMPO_EA_Right(StateComponent& F1, StateComponent const& F0,
                             LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                             LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                             BasicTriangularMPO const& Op, MatrixOperator const& TLeft,
                             MatrixOperator const& TRight, std::complex<double> ExpIK,
                             double Tol = DefaultTol, int Verbose = 0);

void
SolveStringMPO_EA_Left(StateComponent& E1, StateComponent const& E0,
                       LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                       LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                       ProductMPO const& Op, MatrixOperator const& TLeft,
                       MatrixOperator const& TRight, std::complex<double> ExpIK,
                       double Tol = DefaultTol, int Verbose = 0);

void
SolveStringMPO_EA_Right(StateComponent& F1, StateComponent const& F0,
                        LinearWavefunction const& PsiLeft, LinearWavefunction const& PsiRight,
                        LinearWavefunction const& PsiTri, QuantumNumber const& QShift,
                        ProductMPO const& Op, MatrixOperator const& TLeft,
                        MatrixOperator const& TRight, std::complex<double> ExpIK,
                        double Tol = DefaultTol, int Verbose = 0);
#endif
