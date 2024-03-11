// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-dmrg-test.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg2.h"
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
      int NumIter = 4;
      int MinStates = 50;
      int MaxStates = 100000;
      double MixFactor = 0;
      bool TwoSite = false;
      int NumSweeps = 4;
      double MinTrunc = 1;
      int EnvMinStates = 0;
      int EnvMaxStates = 50;
      double EnvTrunc = 1;
      double EnvMixFactor = 0;
      double DiscardMixFactor = 0;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "wavefunction to apply DMRG (required)")
         ("iter,i", prog_opt::value<int>(&NumIter), "Number of Lanczos iterations per step [default 4]")
         ("min-states,m", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 50]")
         ("max-states,x", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 100000]")
         ("trunc,t", prog_opt::value<double>(&MinTrunc),
          ("Minimum permitted truncation error [default "
          +boost::lexical_cast<std::string>(MinTrunc)+"]").c_str())
         ("env-min-states", prog_opt::value(&EnvMinStates),
          "Minimum number of states in the environment [default 0]")
         ("env-max-states", prog_opt::value(&EnvMaxStates),
          "Minimum number of states in the environment [default 50]")
         ("env-weight", prog_opt::value(&EnvTrunc),
          "Kept weight for the environment states [default trunc]")
         ("mix-factor,f", prog_opt::value<double>(&MixFactor),
          "Mixing coefficient for the density matrix [default 0]")
         ("env-mix-factor", prog_opt::value<double>(&EnvMixFactor),
          "Mixing coefficient for the environment density matrix [default 0]")
         ("discard-mix-factor", prog_opt::value<double>(&DiscardMixFactor),
          "Mixing coefficient for the discarded states (large value gives 2-site DMRG) [default 0]")
         ("sweeps,s", prog_opt::value<int>(&NumSweeps), "Number of full sweeps to perform [default 2]")
         //         ("orthogonal", prog_opt::value<std::vector<std::string> >(),
         //          "force the wavefunction to be orthogonal to this state")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-dmrg [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout << "Starting DMRG...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      // Open the wavefunction
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<LinearWavefunction> P = pheap::OpenPersistent(InputWavefunction, CacheSize);
      LinearWavefunction Psi = *P;
      P = pvalue_ptr<LinearWavefunction>();

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
         Psi.Attributes()["Hamiltonian"] = HamString;
      }
      else
      {
         if (Psi.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the initial wavefunction.\n";
            return 1;
         }
         HamString = Psi.Attributes()["Hamiltonian"].as<std::string>();
      }
      std::cout << "Hamiltonian: " << HamString << std::endl;
      MPOperator Hamiltonian = ParseOperator(HamString);

      // Now we can construct the actual DMRG object
      DMRG dmrg(Psi, Hamiltonian);
      Psi = LinearWavefunction(); // we don't need Psi anymore, it will take up space on disk

      // Default value for env-trunc
      if (vm.count("env-trunc") == 0)
         EnvTrunc = MinTrunc;

      dmrg.NumIterations = NumIter;
      dmrg.MinStates = MinStates;
      dmrg.MaxStates = MaxStates;
      dmrg.MinTrunc = MinTrunc;
      dmrg.EnvMinStates = EnvMinStates;
      dmrg.EnvMaxStates = EnvMaxStates;
      dmrg.EnvKeptWeight = EnvTrunc;
      dmrg.MixFactor = MixFactor;
      dmrg.EnvMixFactor = EnvMixFactor;
      dmrg.DiscardMixFactor = DiscardMixFactor;

      double MinTruncPerStep = MinTrunc / (Psi.size()-1);

      int Step = 0;
      for (int Sweep = 0; Sweep < NumSweeps; ++Sweep)
      {
         while (dmrg.Hi != dmrg.Ham.end())
         {
            dmrg.DoIterationMoveRight();
            std::cout << "Step: " << Step++ << '\n';
         }
         dmrg.RightEndpoint();
         while (dmrg.Hi != dmrg.Ham.begin())
         {
            dmrg.DoIterationMoveLeft();
            std::cout << "Step: " << --Step << '\n';
         }
         dmrg.LeftEndpoint();
      }

      //      std::cout << "Finished." << std::endl;
      //      P = pvalue_ptr<LinearWavefunction>(new LinearWavefunction(dmrg.Wavefunction()));
      //     pheap::ShutdownPersistent(P);

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
