// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-random.cpp
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

#include "mp-algorithms/random_wavefunc.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/proccontrol.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "lattice/infinitelattice.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      double Beta = 3;
      int Count = 20;
      std::string LatticeFile;
      std::string FName;
      std::string Target;
      int Size;
      bool Force = false;
      bool Quiet = false;
      int Verbose = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help,h", "show this help message")
         ("lattice,l", prog_opt::value(&LatticeFile), "Lattice file (required)")
         ("unitcell,u", prog_opt::value(&Size), "Number of sites")
         ("quantumnumber,q", prog_opt::value(&Target), "Target quantum number")
         ("count,c", prog_opt::value(&Count), FormatDefault("Count of m=1 states to make a superposition", Count).c_str())
         ("out,o", prog_opt::value(&FName), "Output file (required)")
         ("beta,b", prog_opt::value(&Beta), FormatDefault("Inverse temperature for monte-carlo sampling", Beta).c_str())
         ("seed,s", prog_opt::value<unsigned int>(),
          ("Random seed [range 0.."+boost::lexical_cast<std::string>(RAND_MAX)+"]").c_str())
         ("force,f", prog_opt::bool_switch(&Force), "Allow overwriting output files")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity (can be used more than once)")
	 ("quiet", prog_opt::bool_switch(&Quiet), "Don't print informational messages")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("lattice") == 0 || vm.count("out") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-random [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX) : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

      QuantumNumber Q(Lattice->GetSymmetryList(), Target);
      QuantumNumber Ident(Lattice->GetSymmetryList());

      std::vector<BasisList> BL = ExtractLocalBasis1(Lattice->GetUnitCell());
      std::vector<BasisList> FullBL = BL;
      while (int(FullBL.size()) < Size)
	 std::copy(BL.begin(), BL.end(), std::back_inserter(FullBL));
      if (Size != int(FullBL.size()))
      {
	 std::cout << "mp-idmrg: fatal: the wavefunction size must be a multiple of the unit cell size\n";
	 return 1;
      }

      if (!Quiet)
	 std::cout << "Working..." << std::flush;
      LinearWavefunction Psi = CreateRandomWavefunction(FullBL, Q, Beta, Ident, Count, Verbose + 1 - Quiet);

      if (!Quiet)
	 std::cout << "Done" << std::endl;

      FiniteWavefunctionLeft PsiL = FiniteWavefunctionLeft::Construct(Psi);
      normalize(PsiL);

      MPWavefunction Wavefunction;
      Wavefunction.Wavefunction() = PsiL;
      Wavefunction.SetDefaultAttributes();

      // History log
      Wavefunction.AppendHistoryNote("Random seed " + std::to_string(RandSeed));
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));

      // save wavefunction
      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(P);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      if (e.Why == "File exists")
         std::cerr << "Note: use --force (-f) option to overwrite.\n";
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
