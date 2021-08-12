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
           std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
           StatesInfo SInfo_, int Verbose_);

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

      // Evolve the chain by one time step using single-site TDVP.
      void Evolve();

      // Evolve the current two-site block and move left.
      void IterateLeft2();

      // Evolve the leftmost two-site block in the chain.
      void EvolveLeftmostSite2();

      // Move right and evolve the next two-site block.
      void IterateRight2();

      // Evolve the chain by one time step using two-site TDVP.
      void Evolve2();

      // Calculate the error measures epsilon_1 and epsilon_2.
      void CalculateEps();

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
      StatesInfo SInfo;
      int Verbose;

      int TStep = 0;

      // Cumulative error measures.
      double Eps1SqSum;
      double Eps2SqSum;
      
      // The maximum bond dimension in the chain (2TDVP only).
      int MaxStates;

      // Cumulative truncation error (2TDVP only).
      double TruncErrSum;
};

#endif
