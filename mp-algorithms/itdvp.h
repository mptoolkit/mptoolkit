// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/itdvp.h
//
// Copyright (C) 2021-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#if !defined(MPTOOLKIT_MP_ALGORITHMS_ITDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_ITDVP_H

#include "tdvp.h"
#include "wavefunction/infinitewavefunctionleft.h"

struct iTDVPSettings : TDVPSettings
{
   double GMRESTol = 1e-13;
   int EvolutionSweeps = 2;
   int NEps = 2;
};

class iTDVP : public TDVP
{
   public:
      iTDVP() {}

      iTDVP(InfiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_, iTDVPSettings Settings_);

      // Return the current wavefunction in left-canonical form.
      InfiniteWavefunctionLeft Wavefunction() const;

      // Orthogonalize the leftmost/rightmost site in the unit cell, performing
      // truncation and post-expansion, and updating the left/right block
      // Hamiltonian.
      void OrthogonalizeLeftmostSite();
      void OrthogonalizeRightmostSite();

      // Evolve LambdaR backwards in time.
      void EvolveLambdaRRight(std::complex<double> Tau);
      void EvolveLambdaRLeft(std::complex<double> Tau);

      // Evolve the chain by sweeping left/right until the FidTol or
      // MaxSweeps is reached.
      void EvolveLeft(std::complex<double> Tau);
      void EvolveRight(std::complex<double> Tau);

      // Evolve the chain by one timestep using single-site TDVP.
      void Evolve(bool Expand);

      // Calculate the error measures epsilon_1 and epsilon_2.
      void CalculateEps();

      // Calculate the error measures epsilon_3 to epsilon_NEps.
      void CalculateEpsN(std::deque<StateComponent> X, std::deque<StateComponent> Y);

      // Pre-expand the left/right environment of the current site.
      void ExpandLeft();
      void ExpandRight();

      // Sweep the unit cell from right to left, pre-expanding the left environments of each site.
      void ExpandBondsLeft();
      // Sweep the unit cell from left to right, pre-expanding the right environments of each site.
      void ExpandBondsRight();

      // Update the Hamiltonian if time-dependent, recalculating the left/right
      // environments.
      void UpdateHamiltonian();
      void UpdateHamiltonianLeft(std::complex<double> t, std::complex<double> dt);
      void UpdateHamiltonianRight(std::complex<double> t, std::complex<double> dt);

      QuantumNumber QShift;
      StateComponent BlockHamL, BlockHamR;

      LinearWavefunction PsiOld;
      MatrixOperator LambdaR;
      MatrixOperator LambdaROld;
      MatrixOperator UBoundary, UBoundaryPrev;
      std::deque<StateComponent> HamLOld, HamROld;
      StateComponent CCenter;

      double GMRESTol;
      int EvolutionSweeps;
      int NEps;

      // Error measures epsilon_3^2 to epsilon_NEps^2.
      std::vector<double> EpsNSqSum;

      // Initial energy per unit cell.
      std::complex<double> InitialE;

      // Current energy per unit cell.
      std::complex<double> E;
};

#endif
