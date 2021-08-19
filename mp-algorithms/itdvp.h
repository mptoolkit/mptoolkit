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

#include "common/leftrightstack.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"

class iTDVP
{
   public:
      iTDVP() {}

      iTDVP(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_,
            std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
            double GMRESTol_, StatesInfo SInfo_, int Verbose_);

      // Return the current wavefunction in left-canonical form.
      InfiniteWavefunctionLeft Wavefunction() const;

      // Update the elements of the effective Hamiltonian.
      void UpdateHam();

      // Evolve the wavefunction by one time step.
      void Evolve();

      // Calculate the error measures epsilon_1 and epsilon_2.
      void CalculateEps();

      // Expand the bond  dimension of the iMPS using the projection of H|Psi>
      // onto the subspace of orthogonal two-site variations.
      void ExpandBonds();

      std::deque<StateComponent> PsiAC;
      std::deque<MatrixOperator> PsiC;
      LinearWavefunction PsiAL;
      LinearWavefunction PsiAR;
      QuantumNumber QShift;
      std::deque<StateComponent> HamL;
      std::deque<StateComponent> HamR;
      BasicTriangularMPO Hamiltonian;

      std::complex<double> Timestep;     // The complex timestep in the form -i*dt.
      int MaxIter;
      double ErrTol;
      StatesInfo SInfo;
      double GMRESTol;
      int Verbose;

      int TStep = 0;

      // Cumulative error measures epsilon_L^2 and epsilon_R^2, which quantify
      // how well the left/right isometric AL/AR correspond to the AC and C
      // matrices.
      double EpsLSqSum = 0.0;
      double EpsRSqSum = 0.0;

      // Cumulative error measures epsilon_1^2 and epsilon_2^2, given by the
      // squared Frobenius norms of the projection of H|Psi> onto the subspace
      // of orthogonal 1- and 2-site variations, respectively.
      double Eps1SqSum = 0.0;
      double Eps2SqSum = 0.0;
      
      // The maximum bond dimension in the chain.
      int MaxStates = 1;
};

#endif
