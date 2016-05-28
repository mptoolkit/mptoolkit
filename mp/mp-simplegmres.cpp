// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-simplegmres.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/simplecg.h"
#include "mp-algorithms/gmres.h"
#include "mp/copyright.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mp-algorithms/solver-gmres.h"
#include <iostream>

int main(int argc, char** argv)
{
   if (argc < 7 || argc > 10)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "GMRES Solver\n";
      std::cerr << "usage: mp-simplegmres <lattice> <operator-prefix> <freq> <broad> <psi> <rhs> "
                << "[num-iter] [<maxstates>] [<correction>]\n";
      return 1;
   }

   DEBUG_TRACE(argc);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[5], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string OpPrefix = argv[2];
   std::string Freq = argv[3];
   std::string Broad = argv[4];
   pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(argv[6]);
   int NumIter = 4;
   if (argc >= 8) NumIter = boost::lexical_cast<int>(argv[7]);
   int MaxStates = DefaultMaxStates;
   if (argc >= 9) MaxStates = boost::lexical_cast<int>(argv[8]);
   double CFactor = 0;
   if (argc >= 10) CFactor = boost::lexical_cast<double>(argv[9]);

   std::string Operator = OpPrefix + "_H(" + Freq + "," + Broad + ")";
   std::string Operator2 = OpPrefix + "2_H(" + Freq + "," + Broad + ")";

   SolverGmres solver(*Psi, (*System)[Operator], (*System)[Operator2], *Rhs);

   std::cout.precision(14);

   MPWavefunction xOld = solver.Wavefunction();
   MPWavefunction AxOld = solver.WavefunctionAx();

   double Eta = boost::lexical_cast<double>(Broad);

   bool First = false;
   for (int Sweeps = 0; Sweeps < 3; ++Sweeps)
   {

   {
      solver.ExpandLeft();
      double E = First ? 0.0 : solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m 
		<< ' ' << (-Eta*solver.NormSq()) << '\n';
   }

   // sweep right
   while (solver.RightSize() > 1)
   {
      solver.ShiftRightAndExpand();
      //solver.ExpandRight();
      double E = First ? 0.0 : solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateLeft(MaxStates, CFactor);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m
		<< ' ' << (-Eta*solver.NormSq()) << '\n';
   }
   First = false;

#if 0
   std::cout << "Residual norm at last iteration: "<< solver.ResidualNorm() <<
	     "\nExact residual norm: " << solver.ExactResidualNorm() << '\n';
   std::cout << "difference x: " 
	     << (difference(solver.Wavefunction(), xOld) / norm_2(xOld))
             << "\ndifference Ax: " 
	     << (difference(solver.WavefunctionAx(), AxOld) / norm_2(AxOld))
             << '\n';
   xOld = solver.Wavefunction();
   AxOld = solver.WavefunctionAx();
#endif

   {
      solver.ExpandRight();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m
		<< ' ' << (-Eta*solver.NormSq()) << '\n';
   }

   // sweep left
   while (solver.LeftSize() > 1)
   {
      solver.ShiftLeftAndExpand();
      //solver.ExpandLeft();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, CFactor);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m
		<< ' ' << (-Eta*solver.NormSq()) << '\n';
   }
   double A2 = solver.ExpectationA2();
   std::cout << "Residual norm at last iteration: "<< solver.ResidualNorm() <<
	     "\nExact residual norm: " << solver.ExactResidualNorm(A2) << '\n';
   std::cout << "difference x: " 
	     << (difference(solver.Wavefunction(), xOld) / norm_2(xOld))
             << "\ndifference Ax: " 
	     << (difference(solver.WavefunctionAx(), AxOld) / norm_2(AxOld))
             << '\n';
   xOld = solver.Wavefunction();
   AxOld = solver.WavefunctionAx();

   }

   *Psi.mutate() = solver.Wavefunction();

   pheap::ShutdownPersistent(Psi);
}
