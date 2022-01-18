// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-tdvp.cpp
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
#include "common/environment.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "lattice/infinite-parser.h"
#include "lattice/infinitelattice.h"
#include "mp-algorithms/ibc-tdvp.h"
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
      std::string BetastepStr;
      std::string InitialTimeStr;
      std::string InitialBetaStr;
      int N = 1;
      int SaveEvery = 1;
      int MaxIter = 20;
      double ErrTol = 1e-16;
      double GMRESTol = 1e-13;
      double FidTol = 1e-12;
      double LambdaTol = 1e-12;
      int MinStates = 1;
      int MaxStates = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = 1e-16;
      bool Expand = false;
      double Eps2SqTol = std::numeric_limits<double>::infinity();
      bool TwoSite = false;
      bool Epsilon = false;
      bool UCExpand = false;
      int NExpand = 0;
      int Verbose = 0;
      int OutputDigits = 0;
      std::string CompositionStr = "secondorder";

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
         ("maxiter", prog_opt::value(&MaxIter),
          FormatDefault("Maximum number of Lanczos iterations per step", MaxIter).c_str())
         ("errtol", prog_opt::value(&ErrTol),
          FormatDefault("Error tolerance for the Lanczos evolution", ErrTol).c_str())
         ("gmrestol", prog_opt::value(&GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", GMRESTol).c_str())
         ("min-states", prog_opt::value(&MinStates),
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("states,m", prog_opt::value(&MaxStates),
          FormatDefault("Maximum number of states", MaxStates).c_str())
         ("trunc,r", prog_opt::value(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("expand,x", prog_opt::bool_switch(&Expand), "Use single-site TDVP with bond dimension expansion")
         ("eps2sqtol,e", prog_opt::value(&Eps2SqTol), "Expand the bond dimension in the next step if Eps2SqSum rises above this value [1TDVP only]")
         ("two-site,2", prog_opt::bool_switch(&TwoSite), "Use two-site TDVP")
         ("fidtol,f", prog_opt::value(&FidTol),
          FormatDefault("Tolerance in the boundary fidelity for expanding the window", FidTol).c_str())
         ("lambdatol,l", prog_opt::value(&LambdaTol),
          FormatDefault("Tolerance in the boundary Lambda matrix for expanding the window", LambdaTol).c_str())
         ("uc-expand", prog_opt::bool_switch(&UCExpand), "Expand the window by whole unit cells rather than single sites")
         ("n-expand", prog_opt::value(&NExpand), "Expand the window manually every n timesteps.")
         ("epsilon", prog_opt::bool_switch(&Epsilon), "Calculate the error measures Eps1SqSum and Eps2SqSum")
         ("composition,c", prog_opt::value(&CompositionStr), FormatDefault("Composition scheme", CompositionStr).c_str())
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
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Available compositions:" << std::endl;
         for (auto const& c : Compositions)
            std::cerr << c.first << " : " << c.second.Description << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      std::cout << "Starting IBC TDVP..." << std::endl;
      std::cout << "Hamiltonian: " << HamStr << std::endl;
      std::cout << "Wavefunction: " << InputFile << std::endl;
      std::cout << "Composition: " << CompositionStr << std::endl;

      // Load the composition scheme.
      Composition Comp;
      for (auto const& c : Compositions)
      {
         if (c.first == CompositionStr)
            Comp = c.second;
      }
      if (Comp.Order == 0)
      {
         std::cerr << "fatal: invalid composition" << std::endl;
         return 1;
      }

      // Open the wavefunction.
      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);

      IBCWavefunction Psi = PsiPtr->get<IBCWavefunction>();

      // Output prefix - if it wasn't specified then use the wavefunction attribute, or
      // fallback to the wavefunction name.
      if (OutputPrefix.empty())
         OutputPrefix = PsiPtr->Attributes()["Prefix"].as<std::string>();

      if (OutputPrefix.empty())
         OutputPrefix = InputFile;

      // Get the initial time & beta from the wavefunction attributes.
      std::complex<double> InitialTime = 0.0;
      if (InitialTimeStr.empty())
      {
         std::string T = PsiPtr->Attributes()["Time"].as<std::string>();
         if (!T.empty())
            InitialTime.real(std::stod(T));
      }
      if (InitialBetaStr.empty())
      {
         std::string B = PsiPtr->Attributes()["Beta"].as<std::string>();
         if (!B.empty())
            InitialTime.imag(-std::stod(B));
      }

      // Allow both timestep and betastep.
      std::complex<double> Timestep = 0.0;
      if (!TimestepStr.empty())
         Timestep += ParseNumber(TimestepStr);
      if (!BetastepStr.empty())
         Timestep += std::complex<double>(0.0, -1.0) * ParseNumber(BetastepStr);

      if (OutputDigits == 0)
      {
         OutputDigits = std::max(formatting::digits(Timestep), formatting::digits(InitialTime));
      }

      // Hamiltonian.
      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;

      // Get the Hamiltonian from the attributes, if it wasn't supplied.
      // Get it from EvolutionHamiltonian, if it exists, or from Hamiltonian.
      if (HamStr.empty())
      {
         if (PsiPtr->Attributes().count("EvolutionHamiltonian") == 0)
         {
            if (PsiPtr->Attributes().count("Hamiltonian") == 0)
            {
               std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute EvolutionHamiltonian." << std::endl;
               return 1;
            }
            HamStr = PsiPtr->Attributes()["Hamiltonian"].as<std::string>();
         }
         else
         {
            HamStr = PsiPtr->Attributes()["EvolutionHamiltonian"].as<std::string>();
         }
      }

      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

      std::cout << "Maximum number of Lanczos iterations: " << MaxIter << std::endl;
      std::cout << "Error tolerance for the Lanczos evolution: " << ErrTol << std::endl;

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      if (Expand || Eps2SqTol != std::numeric_limits<double>::infinity() || TwoSite)
         std::cout << SInfo << std::endl;

      // Make sure calculation of Eps2Sq is turned on if we are using it to
      // determine whether to expand the bonds.
      if (Eps2SqTol != std::numeric_limits<double>::infinity())
         Epsilon = true;

      IBC_TDVP tdvp(Psi, HamMPO, std::complex<double>(0.0, -1.0)*Timestep,
                    Comp, MaxIter, ErrTol, SInfo, Epsilon, Verbose, GMRESTol,
                    FidTol, LambdaTol, UCExpand, NExpand);

      if (SaveEvery == 0)
         SaveEvery = N;

      // Calculate initial values of epsilon_1 and epsilon_2.
      if (Epsilon)
         tdvp.CalculateEps();

      std::cout << "Timestep=" << 0
                << " Time=" << formatting::format_digits(InitialTime, OutputDigits)
                << " WindowSize=" << tdvp.Psi.size()
                << " MaxStates=" << tdvp.MaxStates
                << " E=" << std::real(tdvp.Energy());

      if (Epsilon)
         std::cout << " Eps1SqSum=" << tdvp.Eps1SqSum
                   << " Eps2SqSum=" << tdvp.Eps2SqSum;

      std::cout << std::endl;

      for (int tstep = 1; tstep <= N; ++tstep)
      {
         if (TwoSite)
            tdvp.Evolve2();
         else
         {
            if (tdvp.Eps2SqSum > Eps2SqTol)
            {
               if (Verbose > 0)
                  std::cout << "Eps2Sq tolerance reached, expanding bond dimension..." << std::endl;
               tdvp.EvolveExpand();
            }
            else if (Expand)
               tdvp.EvolveExpand();
            else
               tdvp.Evolve();
         }

         std::cout << "Timestep=" << tstep
                   << " Time=" << formatting::format_digits(InitialTime+double(tstep)*Timestep, OutputDigits)
                   << " WindowSize=" << tdvp.Psi.size()
                   << " MaxStates=" << tdvp.MaxStates
                   << " E=" << std::real(tdvp.Energy());

         if (TwoSite)
            std::cout << " TruncErrSum=" << tdvp.TruncErrSum;

         if (Epsilon)
            std::cout << " Eps1SqSum=" << tdvp.Eps1SqSum
                      << " Eps2SqSum=" << tdvp.Eps2SqSum;

         std::cout << std::endl;

         // Save the wavefunction.
         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            MPWavefunction Wavefunction;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), OutputDigits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), OutputDigits);
            IBCWavefunction PsiSave = tdvp.Wavefunction();
            Wavefunction.Wavefunction() = std::move(PsiSave);
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
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
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
