// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/dmrg.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_DMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_DMRG_H

#include "common/leftrightstack.h"
#include "mp-algorithms/eigensolver.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"
#include "mp-algorithms/expansion.h"

// The members that are actually used by the DMRG code:
// Basis1TotalDimension
// Basis2TotalDimension
// EndIteration
// ExpandLeftEnvironment
// ExpandRightEnvironment
// Oversampling
// PostExpansionAlgo
// PreExpansionAlgo
// ProjectTwoSiteTangent
// Site
// Solve
// Solver
// StartIteration
// StartSweep
// TruncateAndShiftLeft
// TruncateAndShiftRight
// UseDGKS
// Wavefunction

class DMRG
{
   public:
      explicit DMRG(int Verbose_);

      // Initialize the state of the system from a left orthogonalized wavefunction.
      // All sites in Psi_ are assumed to be left-orthogonal except for the right-most site.
      void InitializeLeftOrtho(LinearWavefunction Psi_, BasicTriangularMPO const& Ham_, StateComponent const& E, StateComponent const& F);

      // Size of the window
      int size() const { return Psi.size(); }

      int Site() const { return Site_; }

      void StartSweep(bool IncrementSweepNumber = true, double Broad = 0);

      void EndSweep();    // statistics for end of sweep

      void StartIteration();  // prepare statistics for start of iteration
      void EndIteration();    // statistics for end of iteration

      LocalEigensolver& Solver() { return Solver_; }

      // Invoke the local eigensolver
      std::complex<double> Solve();

      // returns the energy, calculated from a matrix-vector product
      std::complex<double> Energy() const;

      // returns the last energy obtained from the eigensolver
      std::complex<double> LastEnergy() const { return Solver_.LastEnergy(); }

      // returns ||psi - psi_prev||^2, where psi_prev is the previous
      // wavefunction (set to psi by StartSweep())
      double FidelityLoss() const;

      // returns the dimension of the Basis1() of the active site
      int Basis1TotalDimension() const { return C->Basis1().total_dimension(); }

      // returns the dimension of the Basis2() of the active site
      int Basis2TotalDimension() const { return C->Basis2().total_dimension(); }

      // get the current wavefunction
      FiniteWavefunctionLeft Wavefunction() const;

      // Coarse-grain the current site with the site on the right, i.e. 2-site DMRG
      std::pair<std::complex<double>, TruncationInfo> SolveCoarseGrainRight(StatesInfo const& SInfo);

      // Coarse-grain the current site with the site on the left, i.e. 2-site DMRG
      std::pair<std::complex<double>, TruncationInfo> SolveCoarseGrainLeft(StatesInfo const& SInfo);

      TruncationInfo TruncateAndShiftLeft(StatesInfo const& SInfo, int ExtraStates, int ExtraStatesPerSector);
      TruncationInfo TruncateAndShiftRight(StatesInfo const& SInfo, int ExtraStates, int ExtraStatesPerSector);

      // The old 3S algorithm
      TruncationInfo TruncateAndShiftLeft3S(StatesInfo const& States, double MixFactor);
      TruncationInfo TruncateAndShiftRight3S(StatesInfo const& States, double MixFactor);

      // Pre-expansion: expand the environment basis to that it contains at least StatesWanted states, if possible.
      // Returns the actual environment size.
      int ExpandLeftEnvironment(int StatesWanted, int ExtraStatesPerSector);
      int ExpandRightEnvironment(int StatesWanted, int ExtraStatesPerSector);

      virtual void check_structure() const;
      virtual void debug_check_structure() const;

      // internal use only
      void ShiftRight(MatrixOperator const& Lambda);
      void ShiftLeft(MatrixOperator const& Lambda);

      left_right_stack<StateComponent> HamMatrices;
      LinearWavefunction Psi;
      LinearWavefunction::iterator C;
      int Site_;                          // Site is the index of the iterator C
      BasicTriangularMPO Hamiltonian;
      BasicTriangularMPO::const_iterator H;
      LocalEigensolver Solver_;

      // The wavefunction from the start of the sweep, projected into the basis of C
      StateComponent SweepC;

      // Basis expansion
      PreExpansionAlgorithm PreExpansionAlgo;
      PostExpansionAlgorithm PostExpansionAlgo;
      OversamplingInfo Oversampling;
      bool ProjectTwoSiteTangent;

      // some global statistics
      int TotalNumSweeps;            // this is incremented every sweep
      int TotalNumIterations;        // total number of times the local eigensolver is called
      int TotalNumMultiplies;        // total number of multiplies used by the local eigensolver

      // Statistics which are valid at the end of the sweep
      double LastSweepTime;          // Elapsed time for the previous sweep
      double LastSweepFidelity;      // Overlap of the wavefunction before and after the previous sweep

      // some statistics for current sweep
      int SweepNumIterations;        // number of iterations in the sweep; normally equal to the number of sites
      int SweepSumStates;            // sum of IterationNumStates
      int SweepMaxStates;            // Maximum number of states across the sweep
      int SweepNumMultiplies;        // number of mat-vec multiplies this sweep
      double SweepMinEnergy;         // lowest energy (real part) seen this sweep
      double SweepTotalTruncation;   // cumulative truncation error this sweep
      double SweepMaxEntropy;        // maximum entropy this sweep
      double SweepStartTime;         // wall time at the start of the sweep
      std::vector<std::complex<double>> SweepEnergyEachIteration;  // energy of each iteration on the current sweep
      double SweepEnergyVariance;    // variance of the energy across the sweep

      // some statistics for current iteration
      int IterationNumMultiplies;
      int IterationNumStates;
      std::complex<double> IterationEnergy;
      double IterationTruncation;
      double IterationEntropy;

      int Verbose;
};

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d);
PStream::ipstream& operator>>(PStream::ipstream& out, DMRG& d);

#endif
