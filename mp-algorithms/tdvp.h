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
#include "lattice/infinitelattice.h"

class Composition
{
   public:
      Composition() : Order(0) {}
      Composition(int Order_, std::string Description_, std::vector<double> Gamma_)
         : Order(Order_), Description(Description_), Gamma(Gamma_)
      {
         CHECK(Gamma.size() % 2 == 0);
         CHECK_CLOSE(std::accumulate(Gamma.begin(), Gamma.end(), 0.0), 1.0);
      }

      int Order;
      std::string Description;
      std::vector<double> Gamma;
};

extern std::map<std::string, Composition> Compositions;

// Class to handle time-dependent Hamiltonian operators.
class Hamiltonian
{
   public:
      Hamiltonian() {}

      // If Size == 0, do not rescale Hamiltonian size.
      Hamiltonian(std::string HamStr, int Size = 0,
                  std::string Magnus = "2", std::string TimeVar = "t");

      // Get then Hamiltonian MPO to evolve from t to t + dt.
      BasicTriangularMPO operator()(std::complex<double> t = 0.0, std::complex<double> dt = 0.0) const;

      void set_size(int Size_);

      bool is_time_dependent() const
      {
         return TimeDependent;
      }

   private:
      InfiniteLattice Lattice;
      std::string HamOperator;
      int Size;
      std::string Magnus;
      std::string TimeVar;
      bool TimeDependent;
      BasicTriangularMPO HamMPO;
};

class TDVP
{
   public:
      TDVP() {}

      TDVP(Hamiltonian const& Ham_, std::complex<double> InitialTime_,
           std::complex<double> Timestep_, Composition Comp_, int MaxIter_,
           double ErrTol_, StatesInfo SInfo_, bool Epsilon_, int Verbose_);

      TDVP(FiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_,
           std::complex<double> InitialTime_, std::complex<double> Timestep_,
           Composition Comp_, int MaxIter_, double ErrTol_, StatesInfo SInfo_,
           bool Epsilon_, int Verbose_);

      // Return the current wavefunction in left-canonical form.
      FiniteWavefunctionLeft Wavefunction() const;

      // Calculate the energy.
      std::complex<double> Energy() const;

      // Evolve the current site.
      void EvolveCurrentSite(std::complex<double> Tau);

      // Move the orthogonality center left/right, evolving the lambda matrix
      // backwards in time.
      void IterateLeft(std::complex<double> Tau);
      void IterateRight(std::complex<double> Tau);

      // Sweep left/right through the chain for timestep Tau.
      void SweepLeft(std::complex<double> Tau);
      void SweepRight(std::complex<double> Tau);

      // Calculate the contribution to epsilon_1/2 for the current site.
      void CalculateEps1();
      void CalculateEps12();

      // The final sweep right, in which we also calculate epsilon_1/2.
      void SweepRightFinal(std::complex<double> Tau);

      // Evolve the chain by one timestep using single-site TDVP.
      void Evolve();

      // Expand the dimension of the left bond of the current site using the
      // projection of H|Psi> onto the subspace of orthogonal two-site
      // variations.
      void ExpandLeftBond();

      // Sweep left through the chain for timestep Tau, expanding the bond dimensions.
      void SweepLeftExpand(std::complex<double> Tau);

      // Evolve the chain by one time step using 1TDVP, expanding the bond
      // dimensions on the first sweep.
      void EvolveExpand();

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

      // Recalculate the left/right Hamiltonian environments: for
      // time-dependent Hamiltonians.
      void RecalculateLeftEnvironment();
      void RecalculateRightEnvironment();

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
      bool Epsilon;
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
