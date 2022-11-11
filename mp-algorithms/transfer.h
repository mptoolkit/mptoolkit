// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/transfer.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// Transfer matrices between different wavefunctions.
// This is difficult to support different size unit cells.  So we do not allow this.  If there is a StringOp
// then this operator is allowed to have a different size unit cell and we take the lcm.
//
// The complications with mismatched wavefunction unit cell sizes is that the quantum number shifts will be different.
// This means that the symmetry sector isn't well defined - if we construct the left eigenvector it will have a different
// symmetry sector to the right eigenvector.  Supporting this is too complicated, so we simply disallow it.
// If required the caller can extend the wavefunctions to the lowest common multiple unit cell size prior to obtaining
// the transfer matrix.
//

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TRANSFER_H)
#define MPTOOLKIT_MP_ALGORITHMS_TRANSFER_H

#include <complex>
#include "wavefunction/linearwavefunction.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "mpo/product_mpo.h"

// Returns the left eigenvector of the transfer matrix,  The wavefunctions must have the same size unit cell.
// The string operator could have a unit cell that divides the wavefunction unit cell.
// The quantum number of StringOp determines the sector to evaluate the eigenvector.
// The returned eigenvector is in Basis1() of Psi1/Psi2.
std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, double tol = 1E-14, int Verbose = 0);

// Version that starts with a specified initial guess vector
std::tuple<std::complex<double>, MatrixOperator>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, MatrixOperator InitialGuess, double tol = 1E-14, int Verbose = 0);

// this version gets the N largest magnitude eigenvectors.  Parameter N comes first to avoid confusion with default
// arguments and get_left_transfer_eigenvector
std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, double tol = 1E-14, int Verbose = 0);

// version that takes an intitial guess
std::tuple<std::vector<std::complex<double>>, std::vector<MatrixOperator>>
get_left_transfer_eigenvectors(int N, LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, MatrixOperator InitialGuess, double tol = 1E-14, int Verbose = 0);

//
// Right eigenvectors
//

// Get the right eigenvector.  The returned eigenvector is in the Basis2() of Psi1/Psi2.
std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, double tol = 1E-14, int Verbose = 0);

std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, ProductMPO const& StringOp, MatrixOperator InitialGuess, double tol = 1E-14, int Verbose = 0);

// Gets the principal left/right eigenpair of the transfer matrix.
// The string operator could have a unit cell that divides the wavefunction unit cell.
// The eigenvectors are normalized such that inner_prod(Left, delta_shift(Right, QShift)) = 1
// The returned left eigenvector is in the Basis1() of Psi1/Psi2.
// The returned right eigenvector is in the Basis2() of Psi1/Psi2.
std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                       ProductMPO const& StringOp,
                       double tol = 1E-14, int Verbose = 0);

// version that uses StringOp == identity in the sector q
std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, QuantumNumber const& q, double tol = 1E-14, int Verbose = 0);

// version that uses StringOp == identity in the identity sector
std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift, double tol = 1E-14, int Verbose = 0);

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(InfiniteWavefunctionLeft const& Psi1, InfiniteWavefunctionLeft const& Psi2,
                       QuantumNumber const& q,
                       double tol = 1E-14, int Verbose = 0);

// get the entire spectrum up to NumEigen eigenvalues.  If LeftVectors or RightVectors is not null, then
// also calculate the left/right eigenvectors.  These are returned in the Basis1() / Basis2() respectively.
LinearAlgebra::Vector<std::complex<double>>
get_spectrum_string(LinearWavefunction const& Psi, QuantumNumber const& QShift,
                    ProductMPO const& StringOp,
                    int NumEigen, double tol = 1e-10,
                    LinearAlgebra::Vector<MatrixOperator>* RightVectors = nullptr,
                    LinearAlgebra::Vector<MatrixOperator>* LeftVectors = nullptr,
                    int ncv = 0, bool Sort = false, int Verbose = 0);

// Match left and right eigenvectors to the correct complex eigenvalues.
// This works by finding a pairing of numerically identical eigenvalues,
// and re-ordering the output arrays so that the corresponding eigenvectors
// are in the same array index.
// If the eigenvalues differ in absolute magnitude by more than tol, then print a warning message to cerr.
void
MatchEigenvectors(int n,
                 LinearAlgebra::Vector<std::complex<double>>& LeftValues,
                 std::vector<std::complex<double>>& LeftVectors,
                 LinearAlgebra::Vector<std::complex<double>>& RightValues,
                 std::vector<std::complex<double>>& RightVectors, double tol = 1e-10, bool IgnoreFinalMismatch = false);

#endif
