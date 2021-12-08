// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/transfer.h
//
// Copyright (C) 2015-2021 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                     ProductMPO const& StringOp,
                     double tol = 1E-14, int Verbose = 0);

// The returned eigenvector is in the Basis2() of Psi1/Psi2.
std::tuple<std::complex<double>, MatrixOperator>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                      ProductMPO const& StringOp,
                      double tol = 1E-14, int Verbose = 0);

// Gets the principal left/right eigenpair of the transfer matrix.
// The string operator could have a unit cell that divides the wavefunction unit cell.
// The eigenvectors are normalized such that inner_prod(Left, delta_shift(Right, QShift)) = 1
// The returned left eigenvector is in the Basis1() of Psi1/Psi2.
// The returned left eigenvector is in the Basis2() of Psi1/Psi2.
std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, QuantumNumber const& QShift,
                       ProductMPO const& StringOp,
                       double tol = 1E-14, int Verbose = 0);

std::tuple<std::complex<double>, MatrixOperator, MatrixOperator>
get_transfer_eigenpair(InfiniteWavefunctionLeft const& Psi1, InfiniteWavefunctionLeft const& Psi2,
                       QuantumNumber const& q,
                       double tol = 1E-14, int Verbose = 0);

#endif
