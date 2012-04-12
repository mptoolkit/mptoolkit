// -*- C++ -*- $Id$

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
