// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/itdvp.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021-2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "tdvp.h"
#include "wavefunction/infinitewavefunctionleft.h"

struct iTDVPSettings : TDVPSettings
{
   double GMRESTol = 1e-13;
   double LambdaTol = 1e-16;
   int MaxSweeps = 10;
   int NEps = 2;
};

class iTDVP : public TDVP
{
   public:
      iTDVP() {}

      iTDVP(InfiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_, iTDVPSettings Settings_);

      // Return the current wavefunction in left-canonical form.
      InfiniteWavefunctionLeft Wavefunction() const;

      // Orthogonalize the leftmost/rightmost site in the unit cell, ensuring
      // that the left and right bases of LambdaR are the same, and calculate
      // the left/right block Hamiltonian.
      void OrthogonalizeLeftmostSite();
      void OrthogonalizeRightmostSite();

      // Evolve LambdaR backwards in time.
      void EvolveLambdaRRight(std::complex<double> Tau);
      void EvolveLambdaRLeft(std::complex<double> Tau);

      // Evolve the chain by sweeping left/right until the LambdaTol or
      // MaxSweeps is reached.
      void EvolveLeft(std::complex<double> Tau);
      void EvolveRight(std::complex<double> Tau);

      // Evolve the chain by one timestep using single-site TDVP.
      void Evolve();

      // Calculate the error measures epsilon_1 and epsilon_2.
      // (If Epsilon == false, this just calculates the X/Y arrays for
      // ExpandBonds.)
      void CalculateEps();

      // Calculate the error measures epsilon_3 to epsilon_NEps.
      void CalculateEpsN();

      // Expand the dimension of the right bond of the current site and move right.
      void ExpandRightBond();

      // Expands the bond dimensions of the chain (must be run after
      // CalculateEps to generate X and Y)
      void ExpandBonds();

      // Update the Hamiltonian if time-dependent, recalculating the left/right
      // environments.
      void UpdateHamiltonianLeft(std::complex<double> t, std::complex<double> dt);
      void UpdateHamiltonianRight(std::complex<double> t, std::complex<double> dt);

      LinearWavefunction PsiOld;
      MatrixOperator LambdaR;
      MatrixOperator LambdaROld;
      QuantumNumber QShift;
      std::deque<StateComponent> HamLOld;
      StateComponent BlockHamL;
      StateComponent BlockHamR;
      LinearWavefunction::iterator COld;

      double GMRESTol;
      double LambdaTol;
      int MaxSweeps;
      int NEps;

      // Partial terms used to calculate epsilon_1/2 and expand the bond
      // dimension of the unit cell.
      std::deque<StateComponent> X;
      std::deque<StateComponent> Y;

      // Flag to determine whether X and Y have been calculated for the current
      // timestep.
      bool XYCalculated = false;

      // Flag for whether the time-dependent Hamiltonian has already been
      // updated, for when the bonds are exanded before performing a left
      // sweep.
      bool HamUpdated = true;

      // Error measures epsilon_3^2 to epsilon_NEps^2.
      std::vector<double> EpsNSqSum;

      // Initial energy per unit cell.
      std::complex<double> InitialE;

      // Current energy per unit cell.
      std::complex<double> E;
};

#endif
