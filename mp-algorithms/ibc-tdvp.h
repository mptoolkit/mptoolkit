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
               double ErrTol_, double GMRESTol_, double LambdaTol_, int NExpand_,
               StatesInfo SInfo_, int Verbose_);

      IBCWavefunction Wavefunction() const;

      // Expand the IBC window by adding one unit cell to the left/right.
      // NB: These functions assume that the window is left/right orthogonal
      // respectively, and that the window does not end in the middle of a unit
      // cell.
      void ExpandWindowLeft();
      void ExpandWindowRight();

      // Calculate the fidelity loss of the left/right edge of the window
      // compared to the semi-infinite boundaries.
      // NB: These functions assume that the window is left/right orthogonal
      // respectively, and that the window does not end in the middle of a unit
      // cell.
      double CalculateFidelityLossLeft();
      double CalculateFidelityLossRight();

      // Calculate the Frobenius norm of the difference of the left/right
      // Lambda matrices of the window against the semi-infinite boundaries.
      double CalculateLambdaDiffLeft();
      double CalculateLambdaDiffRight();

      // Sweep left/right, expanding the window if LambdaDiff exceeds LambdaTol.
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
      double LambdaTol;
      int NExpand;

      InfiniteWavefunctionLeft PsiLeft;
      InfiniteWavefunctionRight PsiRight;
      BasicTriangularMPO HamiltonianLeft;
      BasicTriangularMPO HamiltonianRight;
      std::deque<StateComponent> HamLeftL;
      std::deque<StateComponent> HamRightR;
};

#endif
