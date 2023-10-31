// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ea-extend.cpp
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp-algorithms/ea-dmrg.h"
#include "lattice/infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      int N = 2;
      std::string HamStr;
      std::string Filename;

      EA_DMRGSettings Settings;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr), "Operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value(&Filename), "EA wavefunction to optimize [required]")
         ("num,n", prog_opt::value(&N), FormatDefault("Number of sweeps", N).c_str())
         ("tol", prog_opt::value(&Settings.Tol), FormatDefault("Error tolerance for the local eigensolver", Settings.Tol).c_str())
         ("gmrestol", prog_opt::value(&Settings.GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", Settings.GMRESTol).c_str())
         ("unityepsilon", prog_opt::value(&Settings.UnityEpsilon),
          FormatDefault("Epsilon value for testing eigenvalues near unity", Settings.UnityEpsilon).c_str())
         ("quiet", prog_opt::bool_switch(&Settings.Quiet), "Reduce output")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -w <psi>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      Settings.Verbose = Verbose;

      // Load the wavefunction.
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(Filename.c_str(), mp_pheap::CacheSize());
      EAWavefunction Psi = PsiPtr->get<EAWavefunction>();

      // Load the Hamiltonian.
      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;

      // Get the Hamiltonian from the attributes, if it wasn't supplied.
      if (HamStr.empty())
      {
         if (PsiPtr->Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute Hamiltonian." << std::endl;
            return 1;
         }
         HamStr = PsiPtr->Attributes()["Hamiltonian"].as<std::string>();
      }
      else
         PsiPtr->Attributes()["Hamiltonian"] = HamStr;

      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);
      if (HamMPO.size() < Psi.left().size())
	 HamMPO = repeat(HamMPO, Psi.left().size() / HamMPO.size());

      EA_DMRG dmrg(Psi, HamMPO, Settings);

      dmrg.SolveCurrentSite();
      for (int n = 0; n < N/2; ++n)
         dmrg.SweepLR();

      PsiPtr.mutate()->Wavefunction() = dmrg.Wavefunction();
      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      return 1;
   }
}
