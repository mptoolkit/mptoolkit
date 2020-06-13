// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/random_wavefunc.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// This defines the class WavefunctionDesc, which is a bit mis-named:
// it defines a configuration of local basis states.
// These states can be used as a building block for Monte Carlo.
// The simplest algorithm is constructing a random configuration within
// a particular symmetry sector.  This is done with the
// CreateRandomWavefunction() function.

#if !defined(MPTOOLKIT_MP_ALGORITHMS_RANDOM_WAVEUNC_H)
#define MPTOOLKIT_MP_ALGORITHMS_RANDOM_WAVEUNC_H

#include "wavefunction/linearwavefunction.h"
#include "quantumnumbers/all_symmetries.h"

LinearWavefunction
CreateRandomWavefunction(std::vector<BasisList> const& Basis,
                         QuantumNumber const& q, double Beta);

LinearWavefunction
CreateRandomWavefunction(std::vector<BasisList> const& Basis,
                         QuantumNumber const& q, double Beta,
                         QuantumNumber const& RightBoundary, int NConfig = 20, int Verbose = 0);

#endif
