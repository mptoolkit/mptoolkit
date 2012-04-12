// -*- C++ -*- $Id$

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
