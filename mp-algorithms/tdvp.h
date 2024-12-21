// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/tdvp.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_TDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_TDVP_H

#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"
#include "lattice/infinitelattice.h"
#include "tdvp-compositions.h"
#include "expansion.h"

// Class to handle time-dependent Hamiltonian operators.
class Hamiltonian
{
   public:
      Hamiltonian() {}

      // If Size == 0, do not rescale Hamiltonian size.
      Hamiltonian(std::string HamStr, int Size = 0,
                  std::string Magnus = "2", std::string TimeVar = "t", int Verbose = 0);

      // Get the Hamiltonian MPO to evolve from t to t + dt.
      BasicTriangularMPO operator()(std::complex<double> t = 0.0, std::complex<double> dt = 0.0) const;

      void set_size(int Size_);

      bool is_time_dependent() const { return TimeDependent; }

   protected:
      InfiniteLattice Lattice;
      std::string HamOperator;
      int Size;
      std::string Magnus;
      std::string TimeVar;
      bool TimeDependent;
      BasicTriangularMPO HamMPO;
      int Verbose;
};

struct TDVPSettings
{
   std::complex<double> InitialTime;
   std::complex<double> Timestep;
   Composition Comp;
   int MaxIter = 10;
   double ErrTol = 1e-16;

   StatesInfo SInfo;

   PreExpansionAlgorithm PreExpansionAlgo;
   double PreExpandFactor = 0.1;
   int PreExpandPerSector = 1;

   PostExpansionAlgorithm PostExpansionAlgo;
   double PostExpandFactor = 0.1;
   int PostExpandPerSector = 1;

   bool ProjectTwoSiteTangent = false;
   OversamplingInfo Oversampling;

   bool Epsilon = false;
   bool Normalize = true;
   int Verbose = 0;
};

class TDVP
{
   public:
      TDVP() {}

      TDVP(Hamiltonian const& Ham_, TDVPSettings const& Settings_);

      TDVP(FiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_, TDVPSettings const& Settings_);

      // Return the current wavefunction in left-canonical form.
      FiniteWavefunctionLeft Wavefunction() const;

      // Calculate the energy.
      std::complex<double> Energy() const;

      // Evolve the current site.
      virtual void EvolveCurrentSite(std::complex<double> Tau);

      // Move the orthogonality center left/right, truncating and
      // post-expanding, and then evolve the lambda matrix backwards in time.
      virtual void IterateLeft(std::complex<double> Tau);
      virtual void IterateRight(std::complex<double> Tau);

      // Pre-expand the left/right environment of the current site.
      void ExpandLeft();
      void ExpandRight();

      // Sweep left/right through the chain for timestep Tau.
      void SweepLeft(std::complex<double> Tau, bool Expand = false);
      void SweepRight(std::complex<double> Tau, bool Expand = false);

      // Calculate the contribution to epsilon_1/2 for the current site.
      void CalculateEps1();
      void CalculateEps12();

      // The final sweep right, in which we also calculate epsilon_1/2.
      void SweepRightFinal(std::complex<double> Tau, bool Expand = false);

      // Evolve the chain by one timestep using single-site TDVP.
      void Evolve(bool Expand);

      // Evolve the current site and move left using 2TDVP.
      void IterateLeft2(std::complex<double> Tau);

      // Evolve the leftmost two-site block in the chain using 2TDVP.
      void EvolveLeftmostSite2(std::complex<double> Tau);

      // Move right and evolve the next site using 2TDVP.
      void IterateRight2(std::complex<double> Tau);

      // Sweep left/right through the chain for timestep Tau using 2TDVP.
      void SweepLeft2(std::complex<double> Tau);
      void SweepRight2(std::complex<double> Tau);

      // Evolve the chain by one time step using 2TDVP.
      void Evolve2();

      // Calculate the error measures epsilon_1 and epsilon_2.
      void CalculateEps();

      // Update the Hamiltonian if time-dependent, recalculating the left/right
      // environments.
      void UpdateHamiltonianLeft(std::complex<double> t, std::complex<double> dt);
      void UpdateHamiltonianRight(std::complex<double> t, std::complex<double> dt);

      LinearWavefunction Psi;
      int Site;                          // The index of the iterator C.
      LinearWavefunction::iterator C;
      Hamiltonian Ham;
      BasicTriangularMPO HamMPO;
      BasicTriangularMPO::const_iterator H;
      std::deque<StateComponent> HamL;
      std::deque<StateComponent> HamR;
      int LeftStop;                      // The site indices where we stop iterations.
      int RightStop;

      std::complex<double> InitialTime;
      std::complex<double> Timestep;
      std::complex<double> Time;
      Composition Comp;                  // The symmetric composition scheme used to perform a timestep.
      int MaxIter;
      double ErrTol;

      StatesInfo SInfo;

      PreExpansionAlgorithm PreExpansionAlgo;
      double PreExpandFactor;
      int PreExpandPerSector;

      PostExpansionAlgorithm PostExpansionAlgo;
      double PostExpandFactor;
      int PostExpandPerSector;

      OversamplingInfo Oversampling;
      bool ProjectTwoSiteTangent;

      bool Epsilon;
      bool Normalize; // Only used for iTDVP at the moment.
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

      // The logarithm of the norm of the state.
      double LogAmplitude = 0.0;
};

// Truncate then expand the left environment of CRight by adding extra states to CLeft.
// Assumes CRight is the current orthogonality center.
void
ExpandLeftEnvironment(StateComponent& CLeft, StateComponent& CRight,
                      StateComponent const& E, StateComponent const& F,
                      OperatorComponent const& HLeft, OperatorComponent const& HRight,
                      StatesInfo SInfo, double ExpandFactor, int ExpandMinStates,
                      int ExpandMinPerSector, int Verbose);

// Truncate then expand the right environment of CLeft by adding extra states to CRight.
// Assumes CLeft is the current orthogonality center.
void
ExpandRightEnvironment(StateComponent& CLeft, StateComponent& CRight,
                      StateComponent const& E, StateComponent const& F,
                      OperatorComponent const& HLeft, OperatorComponent const& HRight,
                      StatesInfo SInfo, double ExpandFactor, int ExpandMinStates,
                      int ExpandMinPerSector, int Verbose);

#endif
