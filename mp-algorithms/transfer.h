// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/transfer.h
//
// Copyright (C) 2021 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TRANSFER_H)
#define MPTOOLKIT_MP_ALGORITHMS_TRANSFER_H

#include <complex>
#include "wavefunction/linearwavefunction.h"
#include "mpo/product_mpo.h"

std::tuple<std::complex<double>, int, StateComponent>
get_left_transfer_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
                     LinearWavefunction const& Psi2, QuantumNumber const& QShift2,
                     ProductMPO const& StringOp,
                     double tol = 1E-14, int Verbose = 0);

std::tuple<std::complex<double>, int, StateComponent>
get_right_transfer_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
                      LinearWavefunction const& Psi2, QuantumNumber const& QShift2,
                      ProductMPO const& StringOp,
                      double tol = 1E-14, int Verbose = 0);


#endif
