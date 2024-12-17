// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-dmrg.cpp
//
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "common/environment.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "lattice/infinite-parser.h"
#include "lattice/infinitelattice.h"
#include "mp-algorithms/ibc-dmrg.h"
#include "mp/copyright.h"
#include "mpo/basic_triangular_mpo.h"
#include "parser/number-parser.h"
#include "pheap/pheap.h"
#include "quantumnumbers/all_symmetries.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include <iostream>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string HamStr;
      std::string HamStrWindow;
      std::string InputFile;
      int N = 1;
      bool Expand = true;
      bool TwoSite = false;
      int Verbose = 0;
      std::string CompositionStr = "secondorder";
      std::string Magnus = "2";
      std::string TimeVar = "t";
      std::string PreExpandAlgo = "rsvd";
      std::string PostExpandAlgo = "rsvd";

      int EvolutionWindowLeft = 0;
      int EvolutionWindowRight = 0;

      IBC_DMRGSettings Settings;

      // Truncation defaults
      Settings.SInfo.MinStates = 1;
      Settings.SInfo.TruncationCutoff = 0;
      Settings.SInfo.EigenvalueCutoff = 1e-16;

      // Oversampling defaults
      Settings.Oversampling = OversamplingInfo(10, 1.0, 1);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "Operator to use for the Hamiltonian (wavefunction attribute \"EvolutionHamiltonian\")")
         ("window-operator,W", prog_opt::value(&HamStrWindow),
          "Operator to use for the window Hamiltonian")
         ("wavefunction,w", prog_opt::value(&InputFile), "Input wavefunction (required)")
         ("num-timesteps,n", prog_opt::value(&N), FormatDefault("Number of timesteps to calculate", N).c_str())
         ("maxiter", prog_opt::value(&Settings.MaxIter),
          FormatDefault("Maximum number of Lanczos iterations per step", Settings.MaxIter).c_str())
         ("errtol", prog_opt::value(&Settings.ErrTol),
          FormatDefault("Error tolerance for the Lanczos evolution", Settings.ErrTol).c_str())
         ("gmrestol", prog_opt::value(&Settings.GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", Settings.GMRESTol).c_str())
         ("min-states", prog_opt::value(&Settings.SInfo.MinStates),
          FormatDefault("Minimum number of states to keep", Settings.SInfo.MinStates).c_str())
         ("states,m", prog_opt::value(&Settings.SInfo.MaxStates),
          FormatDefault("Maximum number of states", Settings.SInfo.MaxStates).c_str())
         ("trunc,r", prog_opt::value(&Settings.SInfo.TruncationCutoff),
          FormatDefault("Truncation error cutoff", Settings.SInfo.TruncationCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&Settings.SInfo.EigenvalueCutoff),
          FormatDefault("Eigenvalue cutoff threshold", Settings.SInfo.EigenvalueCutoff).c_str())
         ("pre", prog_opt::value(&PreExpandAlgo), FormatDefault("Pre-expansion algorithm, choices are " + PreExpansionAlgorithm::ListAvailable(), PreExpandAlgo).c_str())
         ("pre-factor", prog_opt::value(&Settings.PreExpandFactor), FormatDefault("Pre-expansion factor", Settings.PreExpandFactor).c_str())
         ("pre-sector", prog_opt::value(&Settings.PreExpandPerSector), "Pre-expansion number of additional states in each quantum number sector [default 1 for fullsvd, rsvd; default 0 for range, random]")
         ("post", prog_opt::value(&PostExpandAlgo), FormatDefault("Post-expansion algorithm, choices are " + PostExpansionAlgorithm::ListAvailable(), PostExpandAlgo).c_str())
         ("post-factor", prog_opt::value(&Settings.PostExpandFactor), FormatDefault("Post-expansion factor", Settings.PostExpandFactor).c_str())
         ("post-sector", prog_opt::value(&Settings.PostExpandPerSector), "Post-expansion number of additional states in each quantum number sector [default 1 for fullsvd, rsvd; default 0 for range, random]")
         ("twositetangent", prog_opt::bool_switch(&Settings.ProjectTwoSiteTangent), "Project onto the two-site tangent space during pre-expansion")
         ("oversample", prog_opt::value(&Settings.Oversampling.Scale), FormatDefault("For random SVD, oversample by this factor", Settings.Oversampling.Scale).c_str())
         ("oversample-min", prog_opt::value(&Settings.Oversampling.Add), FormatDefault("For random SVD, minimum amount of oversampling", Settings.Oversampling.Add).c_str())
         //("two-site,2", prog_opt::bool_switch(&TwoSite), "Use two-site DMRG") // Not implemented yet.
         ("fidtol,f", prog_opt::value(&Settings.FidTol),
          FormatDefault("Tolerance in the boundary fidelity for expanding the window", Settings.FidTol).c_str())
         ("n-expand", prog_opt::value(&Settings.NExpand), "Expand the window manually every n timesteps")
         ("comoving", prog_opt::value(&Settings.Comoving), "Use a comoving window of fixed width of the specified number of sites")
         ("ewleft", prog_opt::value(&EvolutionWindowLeft), "Leftmost site of the initial evolution window (wavefunction attribute \"EvolutionWindowLeft\")")
         ("ewright", prog_opt::value(&EvolutionWindowRight), "Rightmost site of the initial evolution window (wavefunction attribute \"EvolutionWindowRight\")")
         ("epsilon", prog_opt::bool_switch(&Settings.Epsilon), "Calculate the error measures Eps1SqSum and Eps2SqSum")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
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
         std::cerr << "usage: " << basename(argv[0]) << " -w <input-psi> [options]\n";
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      std::cout << "Starting IBC DMRG...\n";
      std::cout << "Hamiltonian: " << HamStr << '\n';
      if (!HamStrWindow.empty())
         std::cout << "Window operator: " << HamStrWindow << '\n';
      std::cout << "Wavefunction: " << InputFile << '\n';

      Settings.Verbose = Verbose;
      Settings.PreExpansionAlgo = PreExpansionAlgorithm(PreExpandAlgo);
      Settings.PostExpansionAlgo = PostExpansionAlgorithm(PostExpandAlgo);

      // Defaults for the pre- and post-expansion per sector
      if (!vm.count("pre-sector"))
      {
         if (Settings.PreExpansionAlgo == PreExpansionAlgorithm::SVD || Settings.PreExpansionAlgo == PreExpansionAlgorithm::RSVD)
         {
            Settings.PreExpandPerSector = 1;
         }
         else if (Settings.PreExpansionAlgo == PreExpansionAlgorithm::RangeFinding || Settings.PreExpansionAlgo == PreExpansionAlgorithm::Random)
         {
            Settings.PreExpandPerSector = 0;
         }
      }
      if (!vm.count("post-sector"))
      {
         if (Settings.PostExpansionAlgo == PostExpansionAlgorithm::SVD || Settings.PostExpansionAlgo == PostExpansionAlgorithm::RSVD)
         {
            Settings.PostExpandPerSector = 1;
         }
         else if (Settings.PostExpansionAlgo == PostExpansionAlgorithm::RangeFinding || Settings.PostExpansionAlgo == PostExpansionAlgorithm::Random)
         {
            Settings.PostExpandPerSector = 0;
         }
      }

      // Print expansion info.
      std::cout << "Pre-expansion algorithm: " << Settings.PreExpansionAlgo.Name();
      if (Settings.PreExpansionAlgo != PreExpansionAlgorithm::NoExpansion)
      {
         std::cout << " with expansion factor: " << Settings.PreExpandFactor
                   << " per sector: " << Settings.PreExpandPerSector;
      }
      std::cout << '\n';
      std::cout << "Post-expansion algorithm: " << Settings.PostExpansionAlgo.Name();
      if (Settings.PostExpansionAlgo != PostExpansionAlgorithm::NoExpansion)
      {
         std::cout << " with expansion factor: " << Settings.PostExpandFactor
                   << " per sector: " << Settings.PostExpandPerSector;
      }
      std::cout << '\n';
      std::cout << Settings.Oversampling;

      // Open the wavefunction.
      //mp_pheap::InitializeTempPHeap();
      //pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(InputFile, mp_pheap::CacheSize());

      IBCWavefunction Psi = PsiPtr->get<IBCWavefunction>();

      // Get the Hamiltonian from the attributes, if it wasn't supplied.
      // Get it from EvolutionHamiltonian, if it exists, or from Hamiltonian.
      if (HamStr.empty())
      {
         if (PsiPtr->Attributes().count("EvolutionHamiltonian") == 0)
         {
            if (PsiPtr->Attributes().count("Hamiltonian") == 0)
            {
               std::cerr << "fatal: No Hamiltonian specified, use -H option or set wavefunction attribute EvolutionHamiltonian." << std::endl;
               return 1;
            }
            HamStr = PsiPtr->Attributes()["Hamiltonian"].as<std::string>();
         }
         else
         {
            HamStr = PsiPtr->Attributes()["EvolutionHamiltonian"].as<std::string>();
         }
      }

      WindowHamiltonian Ham(HamStr, HamStrWindow, Psi.left().size(), Magnus, TimeVar, Verbose);

      // Get EvolutionWindowLeft/Right from the wavefunction attributes if they
      // weren't specified, or else set them to the window boundaries.
      if (vm.count("ewleft") == 0)
      {
         if (PsiPtr->Attributes().count("EvolutionWindowLeft"))
            EvolutionWindowLeft = PsiPtr->Attributes()["EvolutionWindowLeft"].as<int>();
         else
            EvolutionWindowLeft = Psi.window_offset();
      }

      if (vm.count("ewright") == 0)
      {
         if (PsiPtr->Attributes().count("EvolutionWindowRight"))
            EvolutionWindowRight = PsiPtr->Attributes()["EvolutionWindowRight"].as<int>();
         else
            EvolutionWindowRight = Psi.window_size() - 1 + Psi.window_offset();
      }

      Settings.EvolutionWindowLeft = EvolutionWindowLeft;
      Settings.EvolutionWindowRight = EvolutionWindowRight;

      std::cout << "Maximum number of Lanczos iterations: " << Settings.MaxIter << '\n';
      std::cout << "Error tolerance for the Lanczos evolution: " << Settings.ErrTol << '\n';

      // Turn off bond expansion if we specify no pre-expansion.
      if (Settings.PreExpansionAlgo == PreExpansionAlgorithm::NoExpansion)
         Expand = false;

      if (Expand || TwoSite)
         std::cout << Settings.SInfo << '\n';

      IBC_DMRG dmrg(Psi, Ham, Settings);

      // Calculate initial values of epsilon_1 and epsilon_2.
      if (Settings.Epsilon)
         dmrg.CalculateEps();

      std::cout << "Timestep=" << 0
                << " WindowSize=" << dmrg.Psi.size()
                << " EvWindowSize=" << dmrg.RightStop-dmrg.LeftStop+1
                << " MaxStates=" << dmrg.MaxStates
                << " E=" << std::real(dmrg.Energy());
      if (Settings.Epsilon)
         std::cout << " Eps1SqSum=" << dmrg.Eps1SqSum
                   << " Eps2SqSum=" << dmrg.Eps2SqSum;
      std::cout << '\n';

      for (int tstep = 1; tstep <= N; ++tstep)
      {
         if (TwoSite)
            dmrg.Evolve2();
         else
            dmrg.Evolve(Expand);

         std::cout << "Timestep=" << tstep
                   << " WindowSize=" << dmrg.Psi.size()
                   << " EvWindowSize=" << dmrg.RightStop-dmrg.LeftStop+1
                   << " MaxStates=" << dmrg.MaxStates
                   << " E=" << std::real(dmrg.Energy());
         if (TwoSite)
            std::cout << " TruncErrSum=" << dmrg.TruncErrSum;
         if (Settings.Epsilon)
            std::cout << " Eps1SqSum=" << dmrg.Eps1SqSum
                      << " Eps2SqSum=" << dmrg.Eps2SqSum;
         std::cout << '\n';
      }

      IBCWavefunction PsiSave = dmrg.Wavefunction();

      // Stream the boundaries, if the input file does.
      PsiSave.set_left_filename(Psi.get_left_filename());
      PsiSave.set_right_filename(Psi.get_right_filename());

      MPWavefunction Wavefunction;
      Wavefunction.Wavefunction() = std::move(PsiSave);
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
      Wavefunction.SetDefaultAttributes();
      Wavefunction.Attributes()["Hamiltonian"] = HamStr;
      Wavefunction.Attributes()["EvolutionWindowLeft"] = dmrg.LeftStop;
      Wavefunction.Attributes()["EvolutionWindowRight"] = dmrg.RightStop;

      // Save wavefunction.
      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(P);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      pheap::Cleanup();
      return 1;
   }
}
