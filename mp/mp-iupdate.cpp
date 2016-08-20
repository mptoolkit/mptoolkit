// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iupdate.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "wavefunction/mpwavefunction.h"
#include "wavefunction/infinitewavefunction_old.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "pheap/pheaperror.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "tensor/wigner_eckart.h"
#include "tensor/regularize.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      bool Force = false;
      std::string InputFile;
      std::string OutputFile;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("force,f", prog_opt::bool_switch(&Force),
          "overwrite the output file, if it exists")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("inpsi", prog_opt::value(&InputFile), "inpsi")
         ("outpsi", prog_opt::value(&OutputFile), "outpsi")
         ;

      prog_opt::positional_options_description p;
      p.add("inpsi", 1);
      p.add("outpsi", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("inpsi") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <input-psi> [output-psi]\n";
         std::cerr << desc << '\n';
         std::cerr << "Update an infinite wavefunction file to the latest version of the file format.\n";
         return 1;
      }

      // load an old file format
      pheap::SetExpectedPageFileVersion(1);

      pvalue_ptr<InfiniteWavefunctionOld> Psi;

      if (OutputFile.empty())
      {
         // re-use the input file as the output file
         Psi = pheap::OpenPersistent(InputFile, mp_pheap::CacheSize());
      }
      else
      {
         // create a new file for output
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         // and load the input wavefunction
         Psi = pheap::ImportHeap(InputFile);
      }

      pvalue_ptr<InfiniteWavefunctionLeft> PsiNew = new InfiniteWavefunctionLeft(Make(*Psi));

      pheap::ShutdownPersistent(PsiNew);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
   }
   catch (pheap::PHeapVersionMismatch&)
   {
      std::cerr << "mp-iupdate: no action taken, the wavefunction is already updated.\n";
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
