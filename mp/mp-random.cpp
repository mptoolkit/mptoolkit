// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-random.cpp
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

#include "mp-algorithms/random_wavefunc.h"
#include "mp/copyright.h"
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>
#include "common/hash.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      double Beta = 3;
      int Count = 10;
      std::string LatticeFile;
      std::string OutFile;
      std::string Target;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help,h", "show this help message")
         ("lattice,l", prog_opt::value(&LatticeFile), "Lattice file (required)")
         ("quantumnumber,q", prog_opt::value(&Target), "Target quantum number")
         ("count,c", prog_opt::value(&Count), "Count of m=1 states to make a superposition [default 10]")
         ("out,o", prog_opt::value(&OutFile), "Output file (required)")
         ("beta,b", prog_opt::value(&Beta), "Inverse temperature for monte-carlo sampling [default 3]")
         ("seed,s", prog_opt::value<unsigned int>(),
          ("Random seed [range 0.."+boost::lexical_cast<std::string>(RAND_MAX)+"]").c_str())
         ;

      prog_opt::positional_options_description p;
      p.add("lattice", 1);
      p.add("quantumnumber", 1);
      p.add("count", 1);
      p.add("out", 1);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("lattice") == 0 || vm.count("out") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-random [options]\n";
         std::cerr << "alternative usage: mp-random <lattice> <quantumnumber> <count> <out>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX) : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pheap::Initialize(OutFile, 1, PageSize, CacheSize);
      pvalue_ptr<OperatorList> OpList = pheap::ImportHeap(LatticeFile);

      QuantumNumber Q(OpList->GetSymmetryList(), Target);

      std::cout << "Working." << std::flush;
      MPWavefunction Psi = CreateRandomWavefunction(OpList->GetLattice(), Q, Beta, 1);
      while (Count > 1)
      {
         std::cout << "." << std::flush;
         MPWavefunction P2 = CreateRandomWavefunction(OpList->GetLattice(), Q, Beta, 1);
         P2 *= 2.0 * (double(rand()) / RAND_MAX) - 1.0;
         Psi = Psi + P2;
         --Count;
      }
      Psi.normalize();
      std::cout << "done" << std::endl;

      Psi.Attributes()["CmdLine"] = cmdline(argc, argv);
      Psi.Attributes()["RandomSeed"] = RandSeed;

      //TRACE(overlap(Psi, Psi));

      pvalue_ptr<MPWavefunction> Ret = new MPWavefunction(Psi);
      pheap::ShutdownPersistent(Ret);

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
