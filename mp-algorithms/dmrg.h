// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/dmrg.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_DMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_DMRG_H

#include "common/leftrightstack.h"
#include "common/conflist.h"
#include "common/namedenum.h"
#include "mp-algorithms/eigensolver.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_triangular_mpo.h"
#include <boost/shared_ptr.hpp>
#include <fstream>

struct PreExpansionTraits
{
   enum Enum { SVD, RSVD, RangeFinding, RangeFindingProject, Random, NoExpansion };
   static constexpr Enum Default = RangeFinding;
   static constexpr char const* StaticName = "pre-expansion algorithm";
   static constexpr std::array<char const*,5> Names = { "svd", "rsvd", "range", "random", "none" };
};

using PreExpansionAlgorithm = NamedEnumeration<PreExpansionTraits>;

struct PostExpansionTraits
{
   enum Enum { SVD, RangeFinding, FastRangeFinding, Random, NoExpansion };
   static constexpr Enum Default = FastRangeFinding;
   static constexpr char const* StaticName = "post-expansion algorithm";
   static constexpr std::array<char const*, 5> Names = { "svd", "rangefinding", "fastrangefinding", "random", "none"};
};

using PostExpansionAlgorithm = NamedEnumeration<PostExpansionTraits>;

class DMRG
{
   public:
      typedef left_right_stack<StateComponent> OperatorStackType;

      DMRG() {}

      DMRG(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_);

      // Adds x to the 'orthogonal set', that we explicitly orthogonalize the
      // wavefunction against
      void AddOrthogonalState(FiniteWavefunctionLeft x);

      void StartSweep(bool IncrementSweepNumber = true, double Broad = 0);

      void EndSweep();    // statistics for end of sweep

      void StartIteration();  // prepare statistics for start of iteration
      void EndIteration();    // statistics for end of iteration

      //   int LeftSize() const { return Psi.LeftSize(); }
      //   int RightSize() const { return Psi.RightSize(); }

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
      int BasisTotalDimension1() const { return C->Basis1().total_dimension(); }

      // returns the dimension of the Basis2() of the active site
      int BasisTotalDimension2() const { return C->Basis2().total_dimension(); }

      // get the current wavefunction
      FiniteWavefunctionLeft Wavefunction() const;

      void PrepareConvergenceTest();
      bool IsConverged() const;

      TruncationInfo TruncateAndShiftLeft(StatesInfo const& SInfo, int ExtraStates, int ExtraStatesPerSector);
      TruncationInfo TruncateAndShiftRight(StatesInfo const& SInfo, int ExtraStates, int ExtraStatesPerSector);

      // Expand the environment basis to that it contains at least StatesWanted states, if possible.
      // Returns the actual environment size.
      int ExpandLeftEnvironment(int StatesWanted, int ExtraStatesPerSector);
      int ExpandRightEnvironment(int StatesWanted, int ExtraStatesPerSector);

      void debug_check_structure() const;

      OperatorStackType HamMatrices;
      LinearWavefunction Psi;
      int Site;                          // Site is the index of the iterator C
      LinearWavefunction::iterator C;
      BasicTriangularMPO Hamiltonian;
      BasicTriangularMPO::const_iterator H;
      int LeftStop;                      // the site numbers where we stop iterations
      int RightStop;

      //      std::vector<CenterWavefunction> Ortho;              // set of wavefunctions that we want to be
      //      left_right_stack<MatrixOperator> PsiOrthoProjector;  // orthogonal to

      boost::optional<double> LastOverlap;
      bool IsPsiConverged;
      bool IsConvergedValid;
      KeepListType KeepList;

      // some global statistics
      int TotalSweepNumber;          // this is incremented every sweep
      int TotalSweepRecNumber;       // this is incremented conditionally
      int TotalNumIterations;
      int TotalNumMultiplies;

      // some statistics, for current sweep
      int SweepNumIterations;
      int SweepSumStates;            // sum of IterationNumStates
      int SweepMaxStates;
      int SweepNumMultiplies;        // number of mat-vec multiplies this sweep
      double SweepEnergy;            // lowest energy seen this sweep
      double SweepTruncation;        // cumulative truncation error this sweep
      double SweepEntropy;           // maximum entropy this sweep
      double SweepStartTime;         // wall time at the start of the sweep
      double SweepTruncatedEnergy;   // sum of (E_0 - E_truncated) over the sweep
      double SweepEnergyError;       // standard error of the energy at each iteration
      double SweepLastMixFactor;     // the last used mix factor, for the .sweep log file

      PreExpansionAlgorithm PreExpansionAlgo;
      PostExpansionAlgorithm PostExpansionAlgo;
      double RangeFindingOverhead;
      bool ProjectTwoSiteTangent;

      // some statistics, for current iteration
      int IterationNumMultiplies;
      int IterationNumStates;
      std::complex<double> IterationEnergy;
      double IterationTruncation;
      double IterationEntropy;
      std::complex<double> IterationEnergyBeforeTrunc;
      std::vector<std::complex<double>> IterationEnergyVec;

      // Log files are not serialized, but initialized by CreateLogFiles or
      // RestoreLogFiles
      boost::shared_ptr<std::ofstream> EnergyLog, DiagLog, SweepLog, CpuLog, DensityLog;

      // Following items are not persistent, but loaded from the configuration

      // for convergence, we require 2 conditions:
      // Overlap < ConvergenceOverlapTruncationScale * SweepTruncation
      // and OverlapDifference < ConvergenceOverlapDifferenceOverlapScale * Overlap
      double ConvergenceOverlapTruncationScale;
      double ConvergenceOverlapDifferenceOverlapScale;
      double ConvergenceSweepTruncMin;
      double ConvergenceOverlapMin;
      bool NormalizeWavefunction; // should we normalize the wavefunction after each truncation?
      bool MixUseEnvironment;
      bool UseDGKS;  // use DGKS correction in the lanczos for orthogonal states
      bool DoUpdateKeepList; // set to true to use the KeepList to ensure quantum number subspaces are represented
      LocalEigensolver Solver_;
      int Verbose;

      double MixFactor;
      StateComponent PsiPrevC;
};

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d);
PStream::ipstream& operator>>(PStream::ipstream& out, DMRG& d);

#endif
