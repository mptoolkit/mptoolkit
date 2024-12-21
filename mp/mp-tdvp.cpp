// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-tdvp.cpp
//
// Copyright (C) 2021-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "common/environment.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "mp-algorithms/tdvp.h"
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
      std::string InputFile;
      std::string OutputPrefix;
      std::string TimestepStr;
      int N = 1;
      int SaveEvery = 1;
      bool Expand = true;
      bool TwoSite = false;
      int Verbose = 0;
      int OutputDigits = 0;
      std::string CompositionStr = "secondorder";
      std::string Magnus = "2";
      std::string TimeVar = "t";
      std::string PreExpandAlgo = "rsvd";
      std::string PostExpandAlgo = "rsvd";

      TDVPSettings Settings;

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
         ("wavefunction,w", prog_opt::value(&InputFile), "Input wavefunction (required)")
         ("output,o", prog_opt::value(&OutputPrefix), "Prefix for saving output files")
         ("timestep,t", prog_opt::value(&TimestepStr), "Timestep (required)")
         ("num-timesteps,n", prog_opt::value(&N), FormatDefault("Number of timesteps to calculate", N).c_str())
         ("save-timesteps,s", prog_opt::value(&SaveEvery), "Save the wavefunction every s timesteps")
         ("maxiter", prog_opt::value(&Settings.MaxIter),
          FormatDefault("Maximum number of Lanczos iterations per step", Settings.MaxIter).c_str())
         ("errtol", prog_opt::value(&Settings.ErrTol),
          FormatDefault("Error tolerance for the Lanczos evolution", Settings.ErrTol).c_str())
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
         ("two-site,2", prog_opt::bool_switch(&TwoSite), "Use two-site TDVP")
         ("epsilon", prog_opt::bool_switch(&Settings.Epsilon), "Calculate the error measures Eps1SqSum and Eps2SqSum")
         ("composition,c", prog_opt::value(&CompositionStr), FormatDefault("Composition scheme", CompositionStr).c_str())
         ("magnus", prog_opt::value(&Magnus), FormatDefault("For time-dependent Hamiltonians, use this variant of the Magnus expansion", Magnus).c_str())
         ("timevar", prog_opt::value(&TimeVar), FormatDefault("The time variable for time-dependent Hamiltonians", TimeVar).c_str())
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("timestep") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " -w <input-psi> -t <timestep> [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Available compositions:\n";
         for (auto const& c : Compositions)
         {
            std::cerr << c.first << " : ";
            if (Verbose > 0)
               // Print the full composition.
               std::cerr << c.second;
            else
               std::cerr << c.second.Description;
            std::cerr << '\n';
         }
         std::cerr << std::flush;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      std::cout << "Starting TDVP...\n";
      std::cout << "Hamiltonian: " << HamStr << '\n';
      std::cout << "Wavefunction: " << InputFile << '\n';
      std::cout << "Composition: " << CompositionStr << '\n';

      Settings.Verbose = Verbose;

      // Load the composition scheme.
      for (auto const& c : Compositions)
      {
         if (c.first == CompositionStr)
            Settings.Comp = c.second;
      }
      if (Settings.Comp.Order == 0)
      {
         std::cerr << "fatal: Invalid composition." << std::endl;
         return 1;
      }

      if (Magnus != "2" && Magnus != "4")
      {
         std::cerr << "fatal: Invalid Magnus scheme." << std::endl;
         return 1;
      }

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
      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);

      FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();

      // Output prefix - if it wasn't specified then use the wavefunction attribute, or
      // fallback to the wavefunction name.
      if (OutputPrefix.empty())
         OutputPrefix = PsiPtr->Attributes()["Prefix"].as<std::string>();

      if (OutputPrefix.empty())
         OutputPrefix = InputFile;

      // Get the initial time & beta from the wavefunction attributes.
      std::complex<double> InitialTime = 0.0;

      std::string T = PsiPtr->Attributes()["Time"].as<std::string>();
      if (!T.empty())
         InitialTime.real(std::stod(T));

      std::string B = PsiPtr->Attributes()["Beta"].as<std::string>();
      if (!B.empty())
         InitialTime.imag(-std::stod(B));

      // Allow both timestep and betastep.
      std::complex<double> Timestep = 0.0;
      if (!TimestepStr.empty())
         Timestep += ParseNumber(TimestepStr);

      OutputDigits = std::max(formatting::digits(Timestep), formatting::digits(InitialTime));

      Settings.InitialTime = InitialTime;
      Settings.Timestep = Timestep;

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

      Hamiltonian Ham(HamStr, Psi.size(), Magnus, TimeVar, Verbose);

      std::cout << "Maximum number of Lanczos iterations: " << Settings.MaxIter << '\n';
      std::cout << "Error tolerance for the Lanczos evolution: " << Settings.ErrTol << '\n';

      // Turn off bond expansion if we specify no pre-expansion.
      if (Settings.PreExpansionAlgo == PreExpansionAlgorithm::NoExpansion)
         Expand = false;

      if (Expand || TwoSite)
         std::cout << Settings.SInfo << '\n';

      TDVP tdvp(Psi, Ham, Settings);

      if (SaveEvery == 0)
         SaveEvery = N;

      // Calculate initial values of epsilon_1 and epsilon_2.
      if (Settings.Epsilon)
         tdvp.CalculateEps();

      std::cout << "Timestep=" << 0
                << " Time=" << formatting::format_digits(InitialTime, OutputDigits)
                << " MaxStates=" << tdvp.MaxStates
                << " E=" << std::real(tdvp.Energy());
      if (Settings.Epsilon)
         std::cout << " Eps1SqSum=" << tdvp.Eps1SqSum
                   << " Eps2SqSum=" << tdvp.Eps2SqSum;
      std::cout << '\n';

      for (int tstep = 1; tstep <= N; ++tstep)
      {
         if (TwoSite)
            tdvp.Evolve2();
         else
            tdvp.Evolve(Expand);

         std::cout << "Timestep=" << tstep
                   << " Time=" << formatting::format_digits(InitialTime+double(tstep)*Timestep, OutputDigits)
                   << " MaxStates=" << tdvp.MaxStates
                   << " E=" << std::real(tdvp.Energy());
         if (TwoSite)
            std::cout << " TruncErrSum=" << tdvp.TruncErrSum;
         if (Settings.Epsilon)
            std::cout << " Eps1SqSum=" << tdvp.Eps1SqSum
                      << " Eps2SqSum=" << tdvp.Eps2SqSum;
         std::cout << '\n';

         // Save the wavefunction.
         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            MPWavefunction Wavefunction;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), OutputDigits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), OutputDigits);
            FiniteWavefunctionLeft PsiL = tdvp.Wavefunction();
            Wavefunction.Wavefunction() = std::move(PsiL);
            Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
            Wavefunction.SetDefaultAttributes();
            Wavefunction.Attributes()["Time"] = TimeStr;
            Wavefunction.Attributes()["Beta"] = BetaStr;
            Wavefunction.Attributes()["Prefix"] = OutputPrefix;
            Wavefunction.Attributes()["EvolutionHamiltonian"] = HamStr;
            std::string FName = OutputPrefix;
            if (std::real(InitialTime + double(tstep)*Timestep) != 0.0)
            {
               FName += ".t" + TimeStr;
            }
            if (std::imag(InitialTime + double(tstep)*Timestep) != 0.0)
            {
               FName += ".b" + BetaStr;
            }
            *PsiPtr.mutate() = std::move(Wavefunction);
            pheap::ExportHeap(FName, PsiPtr);

            if (Verbose > 0)
               std::cout << "Wavefunction saved" << std::endl;
         }
      }
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
