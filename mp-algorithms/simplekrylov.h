// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/simplekrylov.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(SIMPLEKRYLOV_HVKJH4378767T76RITFEW8REWH)
#define SIMPLEKRYLOV_HVKJH4378767T76RITFEW8REW
#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/splitoperator.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/operatorstack.h"
#include <boost/optional.hpp>
#include <boost/shared_ptr.hpp>

struct SimpleKrylov
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;
   typedef CenterWavefunction WavefunctionType;

   SimpleKrylov() {}

   // Constructor that takes a guess wavefunction
   SimpleKrylov(CenterWavefunction const& Psi_,
                SplitOperator const& Op_,
                std::complex<double> Timestep_,
                CenterWavefunction const& Guess_,
                bool UsePsi1HPsi0_ = false);

   // Constructor where we use Psi0_ itself as the guess, and we already have the
   // superblock matrix elements
   SimpleKrylov(CenterWavefunction const& Psi0_,
                SplitOperator const& Op_,
                std::complex<double> Timestep_,
                SuperblockOperator const& Psi0_H_Psi0,
                bool UsePsi1HPsi0_ = false);

   virtual ~SimpleKrylov();

   void StartSweep(bool IncrementSweepNumber = true);  // prepare statistics for start of sweep
   void EndSweep();    // statistics for end of sweep

   void StartIteration();  // prepare statistics for start of iteration
   void EndIteration();    // statistics for end of iteration

   int LeftSize() const { return Psi1.LeftSize(); }
   int RightSize() const { return Psi1.RightSize(); }

   // Determines the time-evolved wavefunction, using NumIterations krylov vectors
   double Solve(int NumIterations, double ErrBound);

   // Does a truncation and shifts the Center matrix to the right.
   // The new left basis is automatically expanded.
   void ShiftRightAndExpand();

   // Does a truncation and shifts the Center matrix to the left.
   // The new right basis is automatically expanded.
   void ShiftLeftAndExpand();

   // 'Expands' the left basis to cover all states.
   void ExpandLeft();

   // 'Expands' the right basis to cover all states.
   void ExpandRight();

   void ShiftLeft();
   void ShiftRight();

   void SetupLogFiles(std::string const& Prefix, bool Truncate = false);

   void PrepareConvergenceTest() {}

   bool IsConverged() const { return false; }

   void SetWavefunctionAttributes() const;

   CenterWavefunction& Wavefunction();
   CenterWavefunction const& Wavefunction() const;

   TruncationInfo TruncateLeft(StatesInfo const& States, double MixFactor);
   TruncationInfo TruncateRight(StatesInfo const& States, double MixFactor);

   // returns Psi0.Center() in the Psi1 basis.  This is cached
   MatrixOperator const& k0() const;
   void Invalidatek0();

   // Set Psi0 to Psi1
   void AdvanceTime();

   // debug sanity check that the wavefunction and operator basis agree
   void DebugCheckBasis() const;

   CenterWavefunction Psi0;             // the left hand side that we are solving for
   SplitOperator Ham;                 // the fixed operator
   CenterWavefunction Psi1;             // the fixed right hand side

   mutable boost::optional<MatrixOperator> k0Cache;

   bool UsePsi1HPsi0;  // if true, use Psi1_H_Psi0 to get the second Krylov vector

   SuperblockOperator Psi1_H_Psi0;     // matrix elements   | psi1 > H < psi0 |
   SuperblockOperator Psi1_H_Psi1;     // matrix elements   | psi1 > H < psi1 |
   TransformOperator Psi1_Psi0;        // matrix elements   | psi1 > < psi0 |

   QuantumNumber Ident;

   std::complex<double> Timestep;

   MatrixOperator OldPsi1Center;  // Set to Psi1.Center() at the start of Solve()

   boost::shared_ptr<std::ofstream> SweepLog, StepLog;

   // some global statistics
   int TotalSweepNumber;          // this is incremented every sweep

   // some statistics, for current sweep

   int SweepNumIterations;
   int SweepSumStates;            // sum of IterationNumStates
   int SweepMaxStates;
   int SweepNumMultiplies;        // number of mat-vec multiplies this sweep
   double SweepTruncation;        // cumulative truncation error this sweep
   double SweepRealTruncation;    // cumulative real truncation error this sweep
   double SweepEntropy;           // maximum entropy of the sweep
   double SweepStartTime;         // wall time at the start of the sweep
   double SweepStartCPUTime;      // CPU time at the start of the sweep
   double SweepSolverVariance;     // sum of IterationSolverVariance - Pythagoras should apply here
   double SweepStDevBound;        // sum of sqrt(IterationVariance) + sqrt(IterationRealTruncation)
   double SweepOverallVariance;   // sum of IterationOverallVariance

   // some statistics, for current iteration
   int IterationNumMultiplies;
   int IterationNumStates;
   double IterationTruncation;
   double IterationEntropy;
   double IterationSolverVariance; // the change in wavefunction this iteration: ||Psi_orig - Psi_new||^2
   double IterationSolverTolerance; // residual error estimate from the Lanczos exponential
   double IterationRealTruncation; // the actual truncation error of wavefunction (not incorporating mixing)
   double IterationOverallVariance; // change in wavefunction, including truncation
   double IterationStartTime;  // wall time at the start of the iteration
};

PStream::opstream& operator<<(PStream::opstream& out, SimpleKrylov const& s);
PStream::ipstream& operator>>(PStream::ipstream& in, SimpleKrylov& s);

#endif
