// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/itdvp.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_ITDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_ITDVP_H

#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"

class iTDVP
{
   public:
      iTDVP() {}

      iTDVP(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
            std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
            double GMRESTol_, int MaxSweeps_, double LambdaTol_, StatesInfo SInfo_,
            int NEps_, int Verbose_);

      // Return the current wavefunction in left-canonical form.
      InfiniteWavefunctionLeft Wavefunction() const;

      // Evolve the current site and move left.
      void IterateLeft();

      // Evolve the leftmost site in the chain.
      void EvolveLeftmostSite();

      // Move right and evolve the next site.
      void IterateRight();

      // Orthogonalize the leftmost/rightmost site in the unit cell, ensuring
      // that the left and right bases of LambdaR are the same, and calculate
      // the left/right block Hamiltonian.
      void OrthogonalizeLeftmostSite();
      void OrthogonalizeRightmostSite();

      // Evolve LambdaR backwards in time.
      void EvolveLambdaRRight();
      void EvolveLambdaRLeft();

      // Evolve the chain by one time step using single-site TDVP.
      void Evolve();

      // Calculate the error measures epsilon_1 and epsilon_2.
      void CalculateEps();

      // Calculate the error measures epsilon_3 to epsilon_NEps.
      void CalculateEpsN();

      // Expand the dimension of the right bond of the current site and move right.
      void ExpandRightBond();

      // Expands the bond dimensions of the chain (must be run after
      // CalculateEps to generate X and Y)
      void ExpandBonds();

      LinearWavefunction Psi;
      MatrixOperator LambdaR;
      QuantumNumber QShift;
      std::deque<StateComponent> HamL;
      std::deque<StateComponent> HamLOld;
      std::deque<StateComponent> HamR;
      StateComponent BlockHamL;
      StateComponent BlockHamR;
      int Site;                          // The index of the iterator C.
      LinearWavefunction::iterator C;
      LinearWavefunction::iterator COld;
      BasicTriangularMPO Hamiltonian;
      BasicTriangularMPO::const_iterator H;
      int LeftStop;                      // The site indices where we stop iterations,
      int RightStop;

      std::complex<double> Timestep;     // The complex timestep in the form -i*dt.
      int MaxIter;
      double ErrTol;
      StatesInfo SInfo;
      double GMRESTol;
      double LambdaTol;
      int MaxSweeps;
      int NEps;
      int Verbose;

      int TStep = 0;

      // Partial terms used to calculate epsilon_1/2 and expand the bond
      // dimension of the unit cell.
      std::deque<StateComponent> X;
      std::deque<StateComponent> Y;

      // Cumulative error measures epsilon_1^2 and epsilon_2^2, given by the
      // squared Frobenius norms of the projection of H|Psi> onto the subspace
      // of orthogonal 1- and 2-site variations, respectively.
      double Eps1SqSum;
      double Eps2SqSum;

      // Error measures epsilon_3^2 to epsilon_NEps^2.
      std::vector<double> EpsNSqSum;

      // The maximum bond dimension in the chain.
      int MaxStates = 1;

      // Initial energy per unit cell.
      std::complex<double> InitialE;

      // Current energy per unit cell.
      std::complex<double> E;
};

#endif
