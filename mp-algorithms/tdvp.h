// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tdvp.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_TDVP_H

#include "common/leftrightstack.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"

class TDVP
{
   public:
      typedef left_right_stack<StateComponent> OperatorStackType;

      TDVP() {}

      TDVP(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
           std::complex<double> Timestep_, int MaxIter_, double ErrTol_, int Verbose_);

      // Return the current wavefunction.
      FiniteWavefunctionLeft Wavefunction() const;

      // Calculate the energy.
      std::complex<double> Energy() const;

      // Evolve the current site and move left.
      void IterateLeft();

      // Evolve the leftmost site in the chain.
      void EvolveLeftmostSite();

      // Move right and evolve the next site.
      void IterateRight();

      // Evolve the chain by one time step.
      void Evolve();

      OperatorStackType HamMatrices;
      LinearWavefunction Psi;
      int Site;                          // The index of the iterator C.
      LinearWavefunction::iterator C;
      BasicTriangularMPO Hamiltonian;
      BasicTriangularMPO::const_iterator H;
      int LeftStop;                      // The site indices where we stop iterations,
      int RightStop;
      std::complex<double> Timestep;    // The complex timestep in the form -i*dt.
      int MaxIter;
      double ErrTol;
      int Verbose;

      // Cumulative error measures.
      double Eps1SqSum;
      double Eps2SqSum;
};

#endif
