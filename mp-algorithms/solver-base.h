// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/solver-base.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

struct SolverBase
{
   SolverBase();

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

   double SweepTruncation;        // cumulative truncation error this sweep
   double SweepEntropy;           // maximum entropy of the sweep

   // some statistics, for current iteration
   int IterationNumMultiplies;
   int IterationNumStates;
   double IterationTruncation;
   double IterationEntropy;
};

PStream::opstream& operator<<(PStream::opstream& out, SolverBase const& d);
PStream::opstream& operator>>(PStream::ipstream& in, SolverBase& d);
