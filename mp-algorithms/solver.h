// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/solver.h
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

#if !defined(DMRG_LOOP_H_HHVERTHYEHYYOIUDLE89P)
#define DMRG_LOOP_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/mpstate.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/matrixproduct.h"
#include "matrixproduct/operatorstack.h"
#include "common/conflist.h"
#include <fstream>
#include <boost/shared_ptr.hpp>

struct TruncationInfo
{
   int m;
   double trunc;
   double entropy;
   double trunc_half;
};

struct DMRG
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;

   typedef MPWavefunction WavefunctionType;

   DMRG() {}

   DMRG(MPWavefunction const& Psi_, MPOperator const& Ham_);

   void StartSweep(bool IncrementSweepNumber = true);  // prepare statistics for start of sweep
   void EndSweep();    // statistics for end of sweep

   void StartIteration();  // prepare statistics for start of iteration
   void EndIteration();    // statistics for end of iteration

   int LeftSize() const { return Psi.LeftSize(); }
   int RightSize() const { return Psi.RightSize(); }

   // Calculates the groundstate, using Iterations number
   // of Lanczos steps
   double Solve(int Iterations);

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

   // returns the energy, calculated from a matrix-vector product
   double Energy();

   MPWavefunction& Wavefunction() { return Psi; }
   MPWavefunction const& Wavefunction() const { return Psi; }

   void PrepareConvergenceTest() {}
   bool IsConverged() const { return false; }

   TruncationInfo TruncateLeft(int MaxStates = DefaultMaxStates, double CFactor = 0);

   TruncationInfo TruncateRight(int MaxStates = DefaultMaxStates, double CFactor = 0);

   void CreateLogFiles(std::string const& BasePath, ConfList const& Conf);

   void RestoreLogFiles(std::string const& BasePath, ConfList const& Conf);

   SuperblockOperator HamMatrices;
   MPWavefunction Psi;
   MPOperator Ham;

   QuantumNumber Ident;

   MatrixOperator RhoLeft, RhoRight;

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

   // some statistics, for current iteration
   int IterationNumMultiplies;
   int IterationNumStates;
   double IterationEnergy;
   double IterationTruncation;

   // Log files are not serialized, but initialized by CreateLogFiles or
   // RestoreLogFiles
   boost::shared_ptr<std::ofstream> EnergyLog, DiagLog, SweepLog, CpuLog, DensityLog;
};

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d);
PStream::ipstream& operator>>(PStream::ipstream& out, DMRG& d);

#endif
