// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-coarsegrain.cpp
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OutputFile;
      std::string InputFile;
      bool Force = false;
      int Coarsegrain = 2;


      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("force,f", prog_opt::bool_switch(&Force), "allow overwriting the output file, if it already exists")
         ("coarsegrain", prog_opt::value(&Coarsegrain), FormatDefault("coarse-grain N-to-1", Coarsegrain).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <input_psi> <output_psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFile == OutputFile)
         PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFile);
      }

      if (PsiPtr->is<InfiniteWavefunctionLeft>())
      {
         InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

         if (Psi.size() % Coarsegrain != 0)
         {
            std::cerr << "Wavefunction size is not a multiple of the coarsegran size, expending...\n";
            Psi = repeat(Psi, statistics::lcm(Psi.size(), Coarsegrain) / Psi.size());
         }

         Psi = coarse_grain(Psi, Coarsegrain, Verbose);

         PsiPtr.mutate()->Wavefunction() = Psi;
      }
      else if (PsiPtr->is<FiniteWavefunctionLeft>())
      {
         FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();

         LinearWavefunction PsiLinear(Psi.base_begin(), Psi.base_end());

         PsiPtr.mutate()->Wavefunction() = FiniteWavefunctionLeft::Construct(coarse_grain(PsiLinear, Coarsegrain));
      }
      else
      {
         std::cerr << "mp-coarsegrain: fatal: unsupported wavefunction type." << std::endl;
         return 1;
      }

      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
