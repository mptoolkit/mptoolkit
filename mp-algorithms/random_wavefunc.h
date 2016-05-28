// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/random_wavefunc.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

struct WavefunctionDesc
{
   std::vector<int> State;
   std::vector<QuantumNumber> Height;

   // Initialize the WavefunctionDesc from the given lattice,
   // using a valid (but otherwise unspecified) default initial state.
   //   WavefunctionDesc(Lattice const& L);
   WavefunctionDesc(std::vector<BasisList> const& L);

   //   WavefunctionDesc(Lattice const& L, QuantumNumber RightBoundary);
   WavefunctionDesc(std::vector<BasisList> const& L, QuantumNumber RightBoundary);

   // Change the configuration at a particular site into NewState.
   // The NewQuantumNumber specifies the height of the new state
   // (in the sense of a Bratelli diagram).  Returns true if the
   // flip was successful, false if the flip was rejected as invalid.
   // A flip is invalid if it would cause any of the heights on the
   // Bratelli diagram to be invalid.  (SU(2) example: if
   // a flip decreases the height at one point, this might cause a
   // height somewhere to the left to become negative, which is not possible.)
   bool Flip(std::vector<BasisList> const& Basis, int Site, int NewState, 
             QuantumNumber const& NewQuantumNumber);

   QuantumNumber const& TransformsAs() const { return Height.front(); }

   // verifies that the configuration is valid
   void CheckValid(std::vector<BasisList> const& L) const;

};

std::ostream& operator<<(std::ostream& out, WavefunctionDesc const& Config);

WavefunctionDesc
CreateRandomConfiguration(std::vector<BasisList> const& Basis, 
                          QuantumNumber const& q, double Beta);

LinearWavefunction 
CreateRandomWavefunction(std::vector<BasisList> const& Basis, 
                         QuantumNumber const& q, double Beta);

LinearWavefunction 
CreateRandomWavefunction(std::vector<BasisList> const& Basis, 
                         QuantumNumber const& q, double Beta,
			 QuantumNumber const& RightBoundary);


#if 0
LinearWavefunction 
CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, double Beta);

LinearWavefunction 
CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, double Beta,
                         QuantumNumber const& RightBoundary);

LinearWavefunction 
CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, double Beta, int Count);

LinearWavefunction 
CreateRandomWavefunction(Lattice const& L, QuantumNumber const& q, double Beta, 
                         QuantumNumber const& RightBoundary, int Count);

// calculate the amplitude of the given configuration in the wavefunction
std::complex<double>
Amplitude(LinearWavefunction const& Psi, WavefunctionDesc const& Config);

#endif

#endif
