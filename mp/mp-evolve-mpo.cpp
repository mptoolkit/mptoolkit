// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-evolve-mpo.cpp
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

#include "common/environment.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/statistics.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "lattice/infinite-parser.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcell-parser.h"
#include "mp/copyright.h"
#include "mpo/expmpo.h"
#include "parser/number-parser.h"
#include "tensor/tensor_eigen.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string HamStr;
      std::string InputFile;
      std::string OutputPrefix;
      std::string TimestepStr;
      std::string Scheme = "1";
      int N = 1;
      int SaveEvery = 1;
      int Verbose = 0;
      int OutputDigits = 0;
      bool Normalize = false;
      bool NoNormalize = false;
      StatesInfo SInfo;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr), "Operator to use for the Hamiltonian (wavefunction attribute \"EvolutionHamiltonian\")")
         ("wavefunction,w", prog_opt::value(&InputFile), "Input wavefunction (required)")
         ("output,o", prog_opt::value(&OutputPrefix), "Prefix for saving output files")
         ("timestep,t", prog_opt::value(&TimestepStr), "Timestep (required)")
         ("num-timesteps,n", prog_opt::value(&N), FormatDefault("Number of timesteps to calculate", N).c_str())
         ("save-timesteps,s", prog_opt::value(&SaveEvery), "Save the wavefunction every s timesteps")
         ("min-states", prog_opt::value(&SInfo.MinStates), FormatDefault("Minimum number of states to keep", SInfo.MinStates).c_str())
         ("max-states", prog_opt::value<int>(&SInfo.MaxStates), FormatDefault("Maximum number of states to keep", SInfo.MaxStates).c_str())
         ("trunc,r", prog_opt::value<double>(&SInfo.TruncationCutoff), FormatDefault("Truncation error cutoff", SInfo.TruncationCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&SInfo.EigenvalueCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", SInfo.EigenvalueCutoff).c_str())
         ("normalize", prog_opt::bool_switch(&Normalize), "Normalize the wavefunction [default true if timestep is real]")
         ("nonormalize", prog_opt::bool_switch(&NoNormalize), "Don't normalize the wavefunction")
         ("scheme,S", prog_opt::value(&Scheme), FormatDefault("Scheme for the MPO exponential (" + ExpMpoScheme::ListAvailable() + ")", Scheme).c_str())
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
         std::cerr << "usage: " << basename(argv[0]) << " -w <input-psi> -t <timestep> [options]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      std::cout << "Starting MPO time evolution..." << std::endl;
      std::cout << "Hamiltonian: " << HamStr << std::endl;
      std::cout << "Wavefunction: " << InputFile << std::endl;

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

      if (Timestep.imag() == 0.0)
         Normalize = true;
      if (NoNormalize)
         Normalize = false;

      if (Normalize)
         std::cout << "Normalizing wavefunction." << std::endl;
      else
         std::cout << "Not normalizing wavefunction." << std::endl;

      OutputDigits = std::max(formatting::digits(Timestep), formatting::digits(InitialTime));

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

      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;
      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

      if (Psi.size() != HamMPO.size())
      {
         int Size = statistics::lcm(Psi.size(), HamMPO.size());
         HamMPO = repeat(HamMPO, Size / HamMPO.size());
      }

      //ProductMPO EvolutionMPO = FirstOrderEvolutionMPO(HamMPO, std::complex<double>(0.0, -1.0) * Timestep);
      ProductMPO EvolutionMPO = expmpo(std::complex<double>(0.0, -1.0) * Timestep * HamMPO, ExpMpoScheme(Scheme));

      // Make the first and last MPO elements a row and column matrix, respectively.
      EvolutionMPO.front() = project_rows(EvolutionMPO.front(), {0});
      EvolutionMPO.back() = project_columns(EvolutionMPO.back(), {0});

      LinearWavefunction PsiL = LinearWavefunction::FromContainer(Psi.begin(), Psi.end());
      std::complex<double> Amplitude = trace(Psi.lambda_l());

      if (SaveEvery == 0)
         SaveEvery = N;

      std::cout << "Timestep=" << 0
                << " Time=" << formatting::format_digits(InitialTime, OutputDigits)
                << std::endl;

      for (int tstep = 1; tstep <= N; ++tstep)
      {
         // Apply the evolution MPO.
         auto MI = EvolutionMPO.begin();
         for (auto I = PsiL.begin(); I != PsiL.end(); ++I, ++MI)
            *I = aux_tensor_prod(*MI, *I);

         std::cout << "Timestep=" << tstep
                   << " Time=" << formatting::format_digits(InitialTime+double(tstep)*Timestep, OutputDigits)
                   << std::endl;

         MatrixOperator M = left_orthogonalize(PsiL, Verbose);
         Amplitude *= trace(M);
         if (Normalize)
            Amplitude *= 1.0 / std::abs(Amplitude);

         truncate_left_orthogonal(PsiL, SInfo, Verbose);

         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            FiniteWavefunctionLeft PsiSave = FiniteWavefunctionLeft::ConstructFromRightOrthogonal(PsiL, Amplitude, Verbose);
            MPWavefunction Wavefunction;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), OutputDigits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), OutputDigits);
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
