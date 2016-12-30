// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/solver-base.cpp
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

#include "solver-base.h"

SolverBase::SolverBase()
   : TotalSweepNumber(0),
     TotalSweepRecNumber(0),
     TotalNumIterations(0),
     TotalNumMultiplies(0)
{
}

PStream::opstream& operator<<(PStream::opstream& out, SolverBase const& d)
{
   return out << d.TotalSweepNumber
              << d.TotalSweepRecNumber
              << d.TotalNumIterations
              << d.TotalNumMultiplies

              << d.SweepNumIterations
              << d.SweepSumStates
              << d.SweepMaxStates
              << d.SweepNumMultiplies
              << d.SweepTruncation
              << d.SweepEntropy

              << d.IterationNumMultiplies
              << d.IterationNumStates
              << d.IterationTruncation
              << d.IterationEntropy
      ;
}

PStream::opstream& operator>>(PStream::ipstream& in, SolverBase& d)
{
   return in >> d.TotalSweepNumber
             >> d.TotalSweepRecNumber
             >> d.TotalNumIterations
             >> d.TotalNumMultiplies

             >> d.SweepNumIterations
             >> d.SweepSumStates
             >> d.SweepMaxStates
             >> d.SweepNumMultiplies
             >> d.SweepTruncation
             >> d.SweepEntropy

             >> d.IterationNumMultiplies
             >> d.IterationNumStates
             >> d.IterationTruncation
             >> d.IterationEntropy
      ;
}
