// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-dmrg-init.cpp
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

#include "pheap/pheap.h"
#include "quantumnumbers/all_symmetries.h"
#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "mp-algorithms/dmrg.h"
#include "mp-algorithms/stateslist.h"
#include "mp-algorithms/dmrgloop.h"
#include "mp/copyright.h"
#include "common/conflist.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "interface/operator-parser.h"
#include <boost/program_options.hpp>
#include <iostream>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      bool Verbose = false;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "initial wavefunction (required)")
         ("config,c", prog_opt::value<std::string>(), "configuration file (required)")
         ("orthogonal", prog_opt::value<std::vector<std::string> >(),
          "force the wavefunction to be orthogonal to this state")
         ("out,o", prog_opt::value<std::string>(),
          "initial part of filename to use for output files (required)")
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "show extra information during the running of mp-dmrg-init")
         ;
      prog_opt::options_description hidden("Hidden options");
      //      hidden.add_options()
      //         ("input-wavefunction", prog_opt::value<std::string>(),
      //          "input wavefunction (required)")
      //         ("base-filename", prog_opt::value<std::string>(),
      //          "base-filename (required)")
      //         ;

      prog_opt::positional_options_description p;
      //      p.add("input-wavefunction", -1);
      //      p.add("base-filename", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("out") == 0
          || vm.count("config") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-dmrg-init [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (vm.count("wavefunction") + vm.count("target") != 1)
      {
         PANIC("either --wavefunction or --target must be specified (but not both!)");
      }

      std::cout << "Starting DMRG.\n";

      std::string ConfigFile = vm["config"].as<std::string>();
      std::cout << "Configuration file: " << ConfigFile << '\n';
      ConfList Conf(ConfigFile);

      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << '\n';

      std::string BasePathFull = vm["out"].as<std::string>();
      std::cout << "Base filename: " << BasePathFull << '\n';
      std::string BasePath, FileName;
      std::tie(BasePath, FileName) = pheap::SplitPathFile(BasePathFull);
      if (BasePath.empty()) BasePath = "./";

      long PageSize = Conf.GetBytes("PageSize", 64*1024);
      long PageCacheSize = Conf.GetBytes("PageCacheSize", 16*1024*1024);
      int NumPageFiles = Conf.Get("NumPageFiles", 1);
      std::string BinPath = Conf.Get("BinPath", std::string("."));
      if (BinPath[BinPath.size()-1] != '/') BinPath += '/';
      std::cout << "Storing checkpoint files in directory: " << BinPath << '\n';

      MessageLogger::Logger("PHEAPFS").SetThreshold(Conf.Get("PHeapLogLevel", 0));

      // copy the configuration file
      std::string sstr = "cp -f " + ConfigFile + " " + BasePath + FileName + ".conf";
      system(sstr.c_str());

      // Make the params file
      std::ofstream Params((BasePath + FileName + ".params").c_str());
      std::copy(argv, argv+argc, std::ostream_iterator<char*>(Params, " "));
      Params << '\n';
      Params.close();

      pheap::Initialize(BinPath+FileName+".bin", NumPageFiles, PageSize, PageCacheSize);

      // load the wavefunction
      MPWavefunction Psi;
      if (vm.count("wavefunction") == 1)
      {
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(InputWavefunction);
         Psi = *P;
      }
      else
         PANIC("fatal: no wavefunction specified.");
      Psi.Attributes()["CmdLine"] = cmdline(argc, argv);

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
      std::cout << "Hamiltonian: " << HamString << '\n';
      MPOperator Hamiltonian = ParseOperator(HamString);

      // Make sure that the wavefunction and the Hamiltonian have the same local basis
#if 0
      CHECK_EQUAL(Psi.size(), Hamiltonian.size());
      for (int i = 0; i < Psi.size(); ++i)
      {
         CHECK_EQUAL(Psi.LookupLinear(i).SiteBasis(), Hamiltonian.LookupLinear(i).SiteBasis());
      }
#endif

      // Make the actual DMRG object
      DMRG dmrg((CenterWavefunction(Psi)), SplitOperator(Hamiltonian), Verbose);

      // Add the orthogonal states
      if (vm.count("orthogonal") != 0)
      {
         std::vector<std::string> OrthoSet = vm["orthogonal"].as<std::vector<std::string> >();
         for (std::size_t j = 0; j < OrthoSet.size(); ++j)
         {
            if (Verbose)
               std::cerr << "Loading orthogonal wavefunction " << OrthoSet[j] << std::endl;

            pvalue_ptr<MPWavefunction> Ortho = pheap::ImportHeap(OrthoSet[j]);
            dmrg.AddOrthogonalState(CenterWavefunction(*Ortho));
            if (!dmrg.Wavefunction().Attributes()["OrthogonalSet"].empty())
               dmrg.Wavefunction().Attributes()["OrthogonalSet"] += " ";
            dmrg.Wavefunction().Attributes()["OrthogonalSet"] += "--orthogonal "
               + quote_shell(OrthoSet[j]);
         }
      }

      std::string NumStatesStr = Conf.Get("NumStates", "");
      StatesList States(NumStatesStr);

      States.ShowOptions(std::cout);

      if (Verbose)
         std::cout << "Creating DMRGLoop object..." << std::endl;

      dmrg.CreateLogFiles(BasePathFull, Conf);
      pvalue_ptr<DMRGLoop<DMRG> > Solver(new DMRGLoop<DMRG>(dmrg, States));

      if (Verbose)
         std::cout << "done." << std::endl;

      pheap::ShutdownPersistent(Solver);
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
