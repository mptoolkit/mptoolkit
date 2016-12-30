// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-add.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/matrixproduct-sum.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int MinStates = 1;
      int MaxStates = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      std::string OutName;
      bool Verbose = false;
      bool Balanced = false;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("min-states", prog_opt::value<int>(&MinStates),
          "Minimum number of states to keep [default 1]")
         ("max-states,m", prog_opt::value<int>(&MaxStates),
          "Maximum number of states to keep [default 100000]")
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          "Cutoff truncation error per site [default 0]")
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          ("Cutoff threshold for density matrix eigenvalues [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("balanced,b", prog_opt::bool_switch(&Balanced),
          "rotate the wavefunctions inputN for N>1, so that the overlap <input1|inputN> is real and positive")
         ("output,o", prog_opt::value(&OutName),
          "output filename [required]")
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "show additional information")
         ;
      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("input", prog_opt::value<std::vector<std::string> >(), "input wavefunction (required)")
         ;

      prog_opt::positional_options_description p;
      p.add("input", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("output") == 0 || vm.count("input") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-add [options] -o output-file-name input1 input2 ...\n";
         std::cerr << desc << "\n";
         return 1;
      }

      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pheap::Initialize(OutName, 1, PageSize, CacheSize);

      std::vector<std::string> InputWavefunctions
         = vm["input"].as<std::vector<std::string> >();

      std::vector<CenterWavefunction> Psi;
      for (unsigned i = 0; i < InputWavefunctions.size(); ++i)
      {
         pvalue_ptr<MPWavefunction> Next = pheap::ImportHeap(InputWavefunctions[i]);
         CenterWavefunction p(*Next);
         if (Balanced && i > 0)
         {
            double c = std::arg(overlap(Psi.front(), p));
            p *= std::polar(1.0, -c);
         }
         Psi.push_back(p);
      }

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      CenterWavefunction TheSum = mpwavefunction_sum(Psi, SInfo, Verbose);

      pvalue_ptr<MPWavefunction> OutPsi = new MPWavefunction(TheSum.AsLinearWavefunction());
      pheap::ShutdownPersistent(OutPsi);
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
