// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-apply-multiple.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "matrixproduct/centerwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/prodoptimizer.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int MinStates = 1;
      int MaxStates = 10000;
      double MinTrunc = 0;
      double MixFactor = 0.01;
      double EigenCutoff = -1;
      bool TwoSite = false;
      int NumSweeps = 2;
      std::string PsiStr, OutStr;
      std::vector<std::string> OperatorStr, RhsStr;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("two-site,2", "modify 2 neighboring sites at once")
         ("min-states", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 1]")
         ("max-states,m", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 10000]")
         ("trunc,r", prog_opt::value<double>(&MinTrunc),
          "Cutoff for the truncation error per site [default 0]")
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          ("Cutoff threshold for density matrix eigenvalues [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("mix-factor,f", prog_opt::value<double>(&MixFactor),
          "Mixing coefficient for the density matrix [default 0.01]")
         ("sweeps,s", prog_opt::value<int>(&NumSweeps), "Number of half-sweeps to perform [default 2]")
         ("wavefunction,w", prog_opt::value(&PsiStr),
          "wavefunction (required)")
         //         ("out,o", prog_opt::value(&OutStr),
         //          "initial part of filename to use for output files (required)")
         ("operator", prog_opt::value(&OperatorStr), "operator list (can supply multiple of these)")
         ("rhs", prog_opt::value(&RhsStr), "right hand side wavefunction (one per operator)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-apply-opt [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      TwoSite = vm.count("two-site");

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(PsiStr, CacheSize);

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = MinTrunc;
      SInfo.EigenvalueCutoff = EigenCutoff;


      std::vector<SplitOperator> OpList;
      std::vector<CenterWavefunction> RhsList;
      for (unsigned i = 0; i < OperatorStr.size(); ++i)
      {
         MPOperator Op = ParseOperator(OperatorStr[i]);
         OpList.push_back(SplitOperator(Op));
         pvalue_ptr<MPWavefunction> ps = pheap::ImportHeap(RhsStr[i]);
         RhsList.push_back(CenterWavefunction(*ps));
      }

      ProductOptimizer Optimizer((CenterWavefunction(*Psi)), OpList, RhsList);
      Psi = pvalue_ptr<MPWavefunction>();
      OpList.clear();
      RhsList.clear();

      for (int Sweeps = 0; Sweeps < NumSweeps; ++Sweeps)
      {
         if (Optimizer.Wavefunction().LeftSize() < Optimizer.Wavefunction().RightSize())
         {
            // sweep right
            Optimizer.ExpandLeft();
            if (TwoSite) Optimizer.ExpandRight();
            Optimizer.Solve();
            TruncationInfo States = Optimizer.TruncateLeft(SInfo, MixFactor);

            // sweep right
            while (Optimizer.RightSize() > 1)
            {
               Optimizer.ShiftRightAndExpand();
               if (TwoSite) Optimizer.ExpandRight();
               double Norm = Optimizer.Solve();
               TruncationInfo States = Optimizer.TruncateLeft(SInfo, MixFactor);

               std::cout << '(' << Optimizer.LeftSize() << ',' << Optimizer.RightSize()
                         << ") " << Norm << ' ' << States.KeptStates() << '\n';
            }
         }
         else
         {
            Optimizer.ExpandRight();
            if (TwoSite) Optimizer.ExpandLeft();
            Optimizer.Solve();
            Optimizer.TruncateRight(SInfo, MixFactor);

            // sweep left
            while (Optimizer.LeftSize() > 1)
            {
               Optimizer.ShiftLeftAndExpand();
               if (TwoSite) Optimizer.ExpandLeft();
               double Norm = Optimizer.Solve();
               TruncationInfo States = Optimizer.TruncateRight(SInfo, MixFactor);

               std::cout << '(' << Optimizer.LeftSize() << ',' << Optimizer.RightSize()
                         << ") " << Norm << ' ' << States.KeptStates() << '\n';
            }
         }
      }
      Psi = new MPWavefunction(Optimizer.Wavefunction().AsLinearWavefunction());
      pheap::ShutdownPersistent(Psi);

   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
