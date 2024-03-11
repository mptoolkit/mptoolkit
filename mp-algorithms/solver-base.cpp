// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/solver-base.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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
