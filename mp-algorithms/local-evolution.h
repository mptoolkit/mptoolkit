// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/local-evolution.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
// Algorithms for local evolution of a wavefunction.  This is essentially TEBD,
// also applies to transfer matrix DMRG etc.
//

#if !defined(LOCAL_EVOLUTION_H_JH7234678HWCIHUIH7234832)
#define LOCAL_EVOLUTION_H_JH7234678HWCIHUIH7234832

#include "matrixproduct/linearwavefunction.h"
#include "matrixproduct/density.h"

// control of the output of SweepRightEvolve and SweepLeftEvolve is via some global variables

extern bool ShowStates;
extern bool ShowEntropy;
extern bool ShowTruncation;
extern bool ShowSmallestKept;

// assuming the input wavefunction is in normal form, applies a list of evolution operators
// to successive bonds.  For a suzuki-trotter decomposition, half of the evolution operators
// will the the identity operation - this should be negligible loss of efficiency as
// the actual application of the evolution term is O(d^2 m^2), versus the O(d^3 m^3)
// singular value decomposition.
// If ShowInfo is true, then information is written to std::cout, as per the global variables
// above.
void SweepRightEvolve(LinearWavefunction& Psi, std::list<SimpleOperator> const& BondOperators,
                      StatesInfo const& SInfo, bool ShowInfo = false);

// Assuming the wavefunction is in right-most normal form, applies a list of bond operators
// from right to left.  The first bond operator to be applied is BondOperators.back(),
// the bond operator at the start of the list is applied last, to the left-most bond.
// At the end of this function, Psi will be in normal form.
void SweepLeftEvolve(LinearWavefunction& Psi, std::list<SimpleOperator> const& BondOperators,
                     StatesInfo const& SInfo, bool ShowInfo = false);

#endif
