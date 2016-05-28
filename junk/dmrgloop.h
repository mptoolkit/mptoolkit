// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// junk/dmrgloop.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "common/statistics.h"
#include "common/messagelogger.h"
#include "conflist.h"
#include "stateslist.h"

template <typename Algorithm>
void DmrgLoop(Algorithm& solver, ConstList const& Conf, StatesList const& States)
{
   int SweepRecNum = 0;
   int SweepNum = 0;

   int const NumIter = Conf.Get("NumIterations", 10);
   msg_log(2, "InfoLog") << "NumIterations = " << NumIter << '\n';

   int const CFactorAvLen = Conf.Get("Correction::MovingAverageLength", 4);
   msg_log(2, "InfoLog") << "Correction::MovingAverageLength = " << CFactorAvLen << '\n';

   double const CFactorMult = Conf.Get("Correction::Multiplier", 10.0);
   msg_log(2, "InfoLog") << "Correction::Multiplier = " << CFactorMult << '\n';

   while (SweepRecNum < States.NumSweeps())
   {
      int MaxStates = States.NumStates(SweepRecNum);
      double CumulativeTruncation = 0;

      if (States.WaitConverge(SweepRecNum))
         solver.PrepareConvergenceTest();

      // The left-most site
      {
         solver.ExpandLeft();
         double E = solver.Solve(NumIter);
         TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
         std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                   << ") " << E << ' ' << States.m << ' ' << States.trunc << '\n';
         CumulativeTruncation += States.trunc;
      }

      // sweep right
      while (solver.RightSize() > 1)
      {
         solver.ShiftRightAndExpand();
         solver.ExpandRight();
         double E = solver.Solve(NumIter);
         TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
         
         std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                   << ") " << E << ' ' << States.m << ' ' << States.trunc << '\n';
         CumulativeTruncation += States.trunc;
      }

      SweepLog << SweepRecNum << ' ' << SweepNum << ' ';

      if (!States.WaitConverge(SweepRecNum) && !States.SaveState(SweepRecNum) 
          && SweepRecNum != States.NumSweeps()-1)
      {
         ++SweepRecNum;
         MaxStates = States.NumStates(SweepRecNum);
      }
      else if (States.WaitConverge(SweepRecNum))
         solver.PrepareConvergenceTest();

      ++SweepNum;

      // The right-most site
      {
         solver.ExpandRight();
         double E = solver.Solve(NumIter);
         TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);
         std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                   << ") " << E << ' ' << States.m << ' ' << States.trunc << '\n';
      }
      
      // sweep left
      while (solver.LeftSize() > 1)
      {
         solver.ShiftLeftAndExpand();
         solver.ExpandLeft();
         double E = solver.Solve(NumIter);
         TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);

         std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
                   << ") " << E << ' ' << States.m << ' ' << States.trunc << '\n';
      }
      
      if (States.SaveState(SweepRecNum) 
          && (!States.WaitConverge(SweepRecNum) || solver.IsConverged()))
      {
         solver.SaveWavefunction(BasePath + BaseFilename + ".psi." 
                                 + boost::lexical_cast<std::string>(SweepRecNum+1));
      }
      
      if (!States.WaitConverge(SweepRecNum) || solver.IsConverged())
         ++SweepRecNum;

      ++SweepNum;
   }

}
