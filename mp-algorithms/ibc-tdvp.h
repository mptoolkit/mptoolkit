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

struct IBC_TDVPSettings : TDVPSettings
{
   double GMRESTol = 1e-13;
   double FidTol = 1e-12;
   double LambdaTol = 1e-12;
   bool UCExpand = false;
   int NExpand = 0;
   int Comoving = 0;
};

class IBC_TDVP : public TDVP
{
   public:
      IBC_TDVP() = default;

      IBC_TDVP(IBCWavefunction const& Psi_, Hamiltonian const& Ham_, IBC_TDVPSettings const& Settings_);

      IBCWavefunction Wavefunction() const;

      // Expand the IBC window by adding one unit cell to the left/right.
      // (Only works if WindowLeft/RightSites == 0.)
      void ExpandWindowLeftUC();
      void ExpandWindowRightUC();

      // Expand the IBC window by adding a site to the left/right.
      void ExpandWindowLeftSite();
      void ExpandWindowRightSite();

      // Expand the IBC window left/right, adding either a site or a unit cell
      // depending on UCExpand.
      void ExpandWindowLeft();
      void ExpandWindowRight();

      // Calculate the fidelity loss of the left/right edge of the window
      // compared to the semi-infinite boundaries.
      double CalculateFidelityLossLeft();
      double CalculateFidelityLossRight();

      // Calculate the Frobenius norm of the difference of the left/right
      // Lambda matrices of the window against the semi-infinite boundaries.
      double CalculateLambdaDiffLeft();
      double CalculateLambdaDiffRight();

      // Sweep left/right, expanding the window if LambdaDiff exceeds FidTol.
      void SweepLeftEW(std::complex<double> Tau, bool Expand);
      void SweepRightEW(std::complex<double> Tau, bool Expand);
      void SweepRightFinalEW(std::complex<double> Tau, bool Expand);

      // Evolve the window for one time step.
      void Evolve(bool Expand);

      // Evolve the window by one time step using 2TDVP.
      void Evolve2();

      double GMRESTol;
      double FidTol;
      double LambdaTol;
      bool UCExpand;
      int NExpand;
      int Comoving;

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
};

#endif
