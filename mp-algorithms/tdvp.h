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

#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"

class TDVP
{
   public:
      TDVP() {}

      TDVP(BasicTriangularMPO const& Ham_, std::complex<double> Timestep_,
           int MaxIter_, double ErrTol_, StatesInfo SInfo_, int Verbose_);

      TDVP(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
           std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
           StatesInfo SInfo_, int Verbose_);

      // Return the current wavefunction in left-canonical form.
      FiniteWavefunctionLeft Wavefunction() const;

      // Calculate the energy.
      std::complex<double> Energy() const;

      // Evolve the current site and move left.
      void IterateLeft();

      // Evolve the leftmost site in the chain.
      void EvolveLeftmostSite();

      // Move right and evolve the next site.
      void IterateRight();

      // Evolve the chain by one time step using single-site TDVP.
      void Evolve();

      // Expand the dimension of the left bond of the current site using the
      // projection of H|Psi> onto the subspace of orthogonal two-site
      // variations.
      void ExpandLeftBond();

      // Evolve the chain by one time step using 1TDVP, expanding the bond
      // dimensions on the right-to-left sweep.
      void EvolveExpand();

      // Evolve the current site and move left using 2TDVP.
      void IterateLeft2();

      // Evolve the leftmost two-site block in the chain using 2TDVP.
      void EvolveLeftmostSite2();

      // Move right and evolve the next site using 2TDVP.
      void IterateRight2();

      // Evolve the chain by one time step using 2TDVP.
      void Evolve2();

      // Calculate the error measures epsilon_1 and epsilon_2.
      void CalculateEps();

      LinearWavefunction Psi;
      int Site;                          // The index of the iterator C.
      LinearWavefunction::iterator C;
      BasicTriangularMPO Hamiltonian;
      BasicTriangularMPO::const_iterator H;
      std::deque<StateComponent> HamL;
      std::deque<StateComponent> HamR;
      int LeftStop;                      // The site indices where we stop iterations,
      int RightStop;

      std::complex<double> Timestep;     // The complex timestep in the form -i*dt.
      int MaxIter;
      double ErrTol;
      StatesInfo SInfo;
      int Verbose;

      int TStep = 0;

      // Cumulative error measures epsilon_1^2 and epsilon_2^2, given by the
      // squared Frobenius norms of the projection of H|Psi> onto the subspace
      // of orthogonal 1- and 2-site variations, respectively.
      double Eps1SqSum;
      double Eps2SqSum;
      
      // The maximum bond dimension in the chain.
      int MaxStates = 1;

      // Cumulative truncation error (2TDVP only).
      double TruncErrSum;
};

#endif
