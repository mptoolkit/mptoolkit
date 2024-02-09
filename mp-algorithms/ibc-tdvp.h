// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-tdvp.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_IBC_TDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_IBC_TDVP_H

#include "tdvp.h"
#include "wavefunction/ibc.h"

class WindowHamiltonian : public Hamiltonian
{
   public:
      WindowHamiltonian() {}

      // If Size == 0, do not rescale Hamiltonian size.
      WindowHamiltonian(std::string HamStrBackground, std::string HamStrWindow,
                        int Size = 0, std::string Magnus = "2", std::string TimeVar = "t", int Verbose = 0);

      // Get the Hamiltonian MPO to evolve from t to t + dt.
      BasicTriangularMPO operator()(int Left, int Right, std::complex<double> t = 0.0, std::complex<double> dt = 0.0) const;

      bool is_window_empty() const { return WindowEmpty; }
      bool is_window_time_dependent() const { return WindowTimeDependent; }

      bool window_size() const { return WindowMPO.size(); }
      bool window_offset() const { return WindowMPO.offset(); }

   protected:
      std::string WindowOperator;
      UnitCellMPO WindowMPO;
      bool WindowEmpty;
      bool WindowTimeDependent;
      int WindowSize;
      int WindowOffset;
};

struct IBC_TDVPSettings : TDVPSettings
{
   double GMRESTol = 1e-13;
   double FidTol = 1e-12;
   int NExpand = 0;
   int Comoving = 0;
   int EvolutionWindowLeft;
   int EvolutionWindowRight;
};

class IBC_TDVP : public TDVP
{
   public:
      IBC_TDVP() = default;

      IBC_TDVP(IBCWavefunction const& Psi_, WindowHamiltonian const& Ham_, IBC_TDVPSettings const& Settings_);

      IBCWavefunction Wavefunction() const;

      // Expand the IBC window by adding a site to the left/right.
      void ExpandWindowLeft();
      void ExpandWindowRight();

      // Expand the section of the window which is being evolved.
      // NOTE: Only works if the current site is the leftmost or rightmost site
      // of the evolution window respectively.
      void ExpandEvolutionWindowLeft();
      void ExpandEvolutionWindowRight();

      // Calculate the fidelity loss of the left/right edge of the window
      // compared to the semi-infinite boundaries.
      double CalculateFidelityLossLeft();
      double CalculateFidelityLossRight();

      // Sweep left/right, expanding the window if FidelityLoss exceeds FidTol.
      void SweepLeftEW(std::complex<double> Tau, bool Expand);
      void SweepRightEW(std::complex<double> Tau, bool Expand);
      void SweepRightFinalEW(std::complex<double> Tau, bool Expand);

      // Evolve the window for one time step.
      void Evolve(bool Expand);

      // Evolve the window by one time step using 2TDVP.
      void Evolve2();

      double GMRESTol;
      double FidTol;
      int NExpand;
      int Comoving;

      WindowHamiltonian HamWindow;
      const InfiniteWavefunctionLeft PsiLeft;
      const InfiniteWavefunctionRight PsiRight;
      QuantumNumber LeftQShift;
      QuantumNumber RightQShift;
      BasicTriangularMPO HamiltonianWindow;
      BasicTriangularMPO HamiltonianLeft;
      BasicTriangularMPO HamiltonianRight;
      std::deque<StateComponent> HamLeftUC;
      std::deque<StateComponent> HamRightUC;
      int WindowLeftSites;
      int WindowRightSites;
      int Offset;

      BasicTriangularMPO::const_iterator HLeft;
      BasicTriangularMPO::const_iterator HRight;
      InfiniteWavefunctionLeft::const_mps_iterator CLeft;
      InfiniteWavefunctionLeft::const_mps_iterator CRight;
      std::deque<StateComponent>::const_iterator HamLeft;
      std::deque<StateComponent>::const_iterator HamRight;

      StateComponent CRefLeft, CRefLeft2;
      StateComponent CRefRight, CRefRight2;
};

#endif
