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

// Generic base class for DMRG.
// There is a tradeoff between making a very generic 'sweeping' algorithm that can be used
// for many purposes, versus having a code that only does one thing but is more readable.
// A very generic algorithm is likely to end up hard to read, with too many customization points.
// This class is intended to cover standard DMRG and iDMRG, including 2-site, EE, and 3S variants,
// and includes all of the ingredients except for the actual 'Sweepleft' / 'SweepRight' functions.
//
// This is a base class that is used by FiniteDMRG and iDMRG.  It is intended that base classes
// will want to override some of the member functions, although mostly they are not declared
// as virtual (there is no need to do this for functions that are not called from the base class).
// There are some virtual functions that are customization points that are called internally,
// ShiftLeft()/ShiftRight() and ModifyLeftBasis()/ModifyRightBasis().
//
// Sketch of the main loop for a DMRG code:
//
// Initialization:
// DMRG dmrg(Verbose);
// dmrg.InitializeLeftOrtho();                   // Initialize the state from a left ortho wavefunction
//
// Right-moving sweep:
// dmrg.StartSweep()                             // reset statistics for the start of a new sweep
// while (dmrg.Site() < dmrg.size()-1)           // loop over sites
// {
//    dmrg.StartIteration();                     // reset statistics for the start of an iteration
//    [ environment pre-expansion, if required ]
//    dmrg.Solve();                              // invoke the solver
//    dmrg.TruncateAndShiftRight();              // truncate the MPS and move the orthogonality center
//    [ display/store statistics, if required ]
//    dmrg.EndIteration();                       // update sweep statistics from the last iteration
// }
// dmrg.EndSweep();                              // update statistics for the full sweep
// [ display/store statistics for the sweep ]
//
// Note that there is no need to explicitly update dmrg.Site(), since this is taken care of
// by the shift of orthogonality center. For a left-to-right sweep, the site index goes from
//  [0,dmrg.size()-1), starting from site 0, and not including the final site. The final site is
// updated when we do a right-to-left sweep, where the site index starts at dmrg.size()-1 down to
// 1.  A left-moving sweep follows the same outline as above, but the loop test is
// while (dmrg.Site() > 0)
// and after calling the solver we need to truncate and shift left.

#if !defined(MPTOOLKIT_MP_ALGORITHMS_DMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_DMRG_H

#include "common/leftrightstack.h"
#include "mp-algorithms/eigensolver.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"
#include "mp-algorithms/expansion.h"

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

      void StartSweep();

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

   protected:
      // internal use only

      // Shift the orthogonality center, updating the hamiltonian matrices, MPO iterator and
      // MPS iterator. On entry *C must satisfy the appropriate orthogonality constraints.
      virtual void ShiftLeft(MatrixOperator const& Lambda);
      virtual void ShiftRight(MatrixOperator const& Lambda);

      // Modify the left and right basis of *C. This is called internally during
      // environment pre-expansion.

      // Modify the left basis of the orthogonality center by the matrix U
      // This is a customization point, allowing derived classes to add additional state
      // that needs updating when the basis changes. The default version simply does
      // *C = prod(U, *C);
      virtual void ModifyLeftBasis(MatrixOperator const& U);

      // Modify the right basis of the orthogonality center by the matrix U
      // This is a customization point, allowing derived classes to add additional state
      // that needs updating when the basis changes. The default version simply does
      // *C = prod(*C, U);
      virtual void ModifyRightBasis(MatrixOperator const& U);

   protected:
      left_right_stack<StateComponent> HamMatrices;
      LinearWavefunction Psi;
      LinearWavefunction::iterator C;
      int Site_;                          // Site is the index of the iterator C
      BasicTriangularMPO Hamiltonian;
      BasicTriangularMPO::const_iterator H;
      LocalEigensolver Solver_;

   public:
      // Basis expansion
      PreExpansionAlgorithm PreExpansionAlgo;
      PostExpansionAlgorithm PostExpansionAlgo;
      OversamplingInfo Oversampling;

      // set to true if pre-expansion should project onto the 2-site trangent space
      // default: false. Don't change unless you know what you are doing.
      bool ProjectTwoSiteTangent;

      // some global statistics. These are updated internally by the {Start|End}{Iteration|Sweep}
      // functions and should otherwise be regarded as read-only.

      int TotalNumSweeps;            // this is incremented every sweep
      int TotalNumIterations;        // total number of times the local eigensolver is called
      int TotalNumMultiplies;        // total number of multiplies used by the local eigensolver

      // Statistics which are valid at the end of the sweep
      double LastSweepTime;          // Elapsed time for the previous sweep

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
