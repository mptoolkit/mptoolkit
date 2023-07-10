// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-itdvp.cpp
//
// Copyright (C) 2004-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "mp-algorithms/itdvp.h"
#include "mp/copyright.h"
#include "mpo/basic_triangular_mpo.h"
#include "parser/number-parser.h"
#include "pheap/pheap.h"
#include "quantumnumbers/all_symmetries.h"
#include "wavefunction/infinitewavefunctionleft.h"
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
      bool Expand = false;
      int Verbose = 0;
      int OutputDigits = 0;
      std::string CompositionStr = "secondorder";
      std::string Magnus = "2";
      std::string TimeVar = "t";
      bool Normalize = false;
      bool NoNormalize = false;

      iTDVPSettings Settings;
      Settings.SInfo.MinStates = 1;
      Settings.SInfo.TruncationCutoff = 0;
      Settings.SInfo.EigenvalueCutoff = 1e-16;

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
         ("lambdatol,l", prog_opt::value(&Settings.LambdaTol),
          FormatDefault("Tolerance for the squared Frobenius norm of the difference of LambdaR for succesive sweeps", Settings.LambdaTol).c_str())
         ("max-sweeps", prog_opt::value(&Settings.MaxSweeps),
          FormatDefault("Maximum number of sweeps", Settings.MaxSweeps).c_str())
         ("epsilon", prog_opt::bool_switch(&Settings.Epsilon), "Calculate the error measures Eps1SqSum and Eps2SqSum")
         ("neps,N", prog_opt::value(&Settings.NEps), "Calculate EpsNSqSum up to N = NEps >= 3")
         ("force-expand", prog_opt::bool_switch(&Settings.ForceExpand), "Force bond dimension expansion [use with caution!]")
         ("composition,c", prog_opt::value(&CompositionStr), FormatDefault("Composition scheme", CompositionStr).c_str())
         ("magnus", prog_opt::value(&Magnus), FormatDefault("For time-dependent Hamiltonians, use this variant of the Magnus expansion", Magnus).c_str())
         ("timevar", prog_opt::value(&TimeVar), FormatDefault("The time variable for time-dependent Hamiltonians", TimeVar).c_str())
         ("normalize", prog_opt::bool_switch(&Normalize), "Normalize the wavefunction [default true if timestep is real]")
         ("nonormalize", prog_opt::bool_switch(&NoNormalize), "Don't normalize the wavefunction")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      // Bond expansion is turned on automatically if an option which uses it
      // is specified, so this option is redundant now. It is kept as a hidden
      // option to ensure compatibility with old scripts.
      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("expand,x", prog_opt::bool_switch(&Expand), "Use single-site TDVP with bond dimension expansion")
         ;

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("timestep") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " -w <input-psi> -t <timestep> [options]" << std::endl;
         std::cerr << desc << std::endl;
         std::cerr << "Available compositions:" << std::endl;
         for (auto const& c : Compositions)
            std::cerr << c.first << " : " << c.second.Description << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      std::cout << "Starting iTDVP..." << std::endl;
      std::cout << "Hamiltonian: " << HamStr << std::endl;
      std::cout << "Wavefunction: " << InputFile << std::endl;
      std::cout << "Composition: " << CompositionStr << std::endl;

      Settings.Verbose = Verbose;

      // Load the composition scheme.
      for (auto const& c : Compositions)
      {
         if (c.first == CompositionStr)
            Settings.Comp = c.second;
      }
      if (Settings.Comp.Order == 0)
      {
         std::cerr << "fatal: Invalid composition" << std::endl;
         return 1;
      }

      // Open the wavefunction.
      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

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

      if (Timestep.imag() == 0.0)
         Normalize = true;
      if (NoNormalize)
         Normalize = false;
      Settings.Normalize = Normalize;

      if (Normalize)
         std::cout << "Normalizing wavefunction." << std::endl;
      else
         std::cout << "Not normalizing wavefunction." << std::endl;

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

      Hamiltonian Ham(HamStr, Psi.size(), Magnus, TimeVar);

      std::cout << "Maximum number of Lanczos iterations: " << Settings.MaxIter << std::endl;
      std::cout << "Error tolerance for the Lanczos evolution: " << Settings.ErrTol << std::endl;

      std::cout << "Maximum number of sweeps: " << Settings.MaxSweeps << std::endl;
      std::cout << "Error tolerance for LambdaR: " << Settings.LambdaTol << std::endl;

      // Turn on bond expansion if trunc or eigen-cutoff have been specified,
      // or if forced bond expansion is specified.
      if (vm.count("trunc") || vm.count("eigen-cutoff") || Settings.ForceExpand)
         Expand = true;

      if (Expand)
         std::cout << Settings.SInfo << std::endl;

      // Make sure calculation of Eps2Sq is turned on if we are calculating EpsNSq for N > 2.
      if (Settings.NEps > 2)
         Settings.Epsilon = true;

      iTDVP itdvp(Psi, Ham, Settings);

      if (SaveEvery == 0)
         SaveEvery = N;

      // Calculate initial values of epsilon.
      if (Settings.Epsilon)
         itdvp.CalculateEps();

      std::cout << "Timestep=" << 0
                << " Time=" << formatting::format_digits(InitialTime, OutputDigits)
                << " MaxStates=" << itdvp.MaxStates
                << " E=" << std::real(itdvp.InitialE);
      if (Settings.Epsilon)
      {
         std::cout << " Eps1SqSum=" << itdvp.Eps1SqSum
                   << " Eps2SqSum=" << itdvp.Eps2SqSum;
         for (int i = 0; i < Settings.NEps-2; ++i)
            std::cout << " Eps" << i+3 << "SqSum=" << itdvp.EpsNSqSum[i];
      }
      std::cout << std::endl;

      for (int tstep = 1; tstep <= N; ++tstep)
      {
         if (Expand)
            itdvp.ExpandBonds();

         itdvp.Evolve();

         std::cout << "Timestep=" << tstep
                   << " Time=" << formatting::format_digits(InitialTime+double(tstep)*Timestep, OutputDigits)
                   << " MaxStates=" << itdvp.MaxStates
                   << " E=" << std::real(itdvp.E);
         if (Settings.Epsilon)
         {
            std::cout << " Eps1SqSum=" << itdvp.Eps1SqSum
                      << " Eps2SqSum=" << itdvp.Eps2SqSum;
            for (int i = 0; i < Settings.NEps-2; ++i)
               std::cout << " Eps" << i+3 << "SqSum=" << itdvp.EpsNSqSum[i];
         }
         std::cout << std::endl;

         // Save the wavefunction.
         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            MPWavefunction Wavefunction;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), OutputDigits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), OutputDigits);
            InfiniteWavefunctionLeft PsiL = itdvp.Wavefunction();
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
