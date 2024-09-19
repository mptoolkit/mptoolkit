// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-normalize.cpp
//
// Copyright (C) 2015-2024 Ian McCulloch <ian@qusim.net>
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
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcell-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      std::string InputFile;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ;

      prog_opt::positional_options_description p;
      p.add("psi1", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi1") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " <psi>\n";
         std::cerr << desc << '\n';
         std::cerr << "Normalizes a finite or IBC wavefunction to 1.0.\n";

         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<MPWavefunction> PsiPtr;
      PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());

      if (PsiPtr->is<FiniteWavefunctionLeft>())
      {
         FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();
         normalize(Psi);
         PsiPtr.mutate()->Wavefunction() = Psi;
      }
      else if (PsiPtr->is<InfiniteWavefunctionLeft>())
      {
         PsiPtr.mutate()->get<InfiniteWavefunctionLeft>().normalize();
      }
      else if (PsiPtr->is<InfiniteWavefunctionRight>())
      {
         PsiPtr.mutate()->get<InfiniteWavefunctionRight>().normalize();
      }
      else
      {
         std::cerr << "mp-normalize: fatal: wavefunction type " << PsiPtr->Type() << " is not supported.\n";
         return 1;
      }

      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
