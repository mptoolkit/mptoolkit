// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-tdvp.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_IBC_TDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_IBC_TDVP_H

#include "tdvp.h"
#include "wavefunction/ibc.h"

class IBC_TDVP : public TDVP
{
   public:
      IBC_TDVP() {}

      IBC_TDVP(IBCWavefunction const& Psi_, BasicTriangularMPO const& Ham_,
               std::complex<double> Timestep_, Composition Comp_, int MaxIter_,
               double ErrTol_, StatesInfo SInfo_, bool Epsilon_, int Verbose_,
               double GMRESTol_, double FidTol_, double LambdaTol_,
               bool UCExpand_ = false, int NExpand_ = 0, int Comoving_ = 0);

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
      void SweepLeftEW(std::complex<double> Tau);
      void SweepRightEW(std::complex<double> Tau);
      void SweepRightFinalEW(std::complex<double> Tau);
      void SweepLeftExpandEW(std::complex<double> Tau);

      // Evolve the window for one time step.
      void Evolve();

      // Evolve the window for one time step, with bond expansion.
      void EvolveExpand();

      // Evolve the window by one time step using 2TDVP.
      void Evolve2();

      double GMRESTol;
      double FidTol;
      double LambdaTol;
      bool UCExpand;
      int NExpand;
      int Comoving;

      InfiniteWavefunctionLeft PsiLeft;
      InfiniteWavefunctionRight PsiRight;
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
