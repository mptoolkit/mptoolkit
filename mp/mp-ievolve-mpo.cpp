// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ievolve-mpo.cpp
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
#include "mpo/basic_triangular_mpo.h"
#include "parser/number-parser.h"
#include "tensor/tensor_eigen.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

ProductMPO
FirstOrderEvolutionMPO(BasicTriangularMPO& HamMPO, std::complex<double> Tau)
{
   std::deque<OperatorComponent> Result;

   for (auto const& I : HamMPO)
   {
      BasisList Basis1(I.Basis1().begin(), I.Basis1().end()-1);
      BasisList Basis2(I.Basis2().begin(), I.Basis2().end()-1);
      OperatorComponent W(I.LocalBasis1(), I.LocalBasis2(), Basis1, Basis2);

      for (int i = 0; i < Basis1.size(); ++i)
      {
         W(i, 0) = I(i, 0) + Tau * I(i, I.size2()-1);

         for (int j = 1; j < Basis2.size(); ++j)
            W(i, j) = I(i, j);
      }

      Result.push_back(W);
   }

   return ProductMPO(GenericMPO(Result.begin(), Result.end()));
}

// Calculate the optimal first-order evolution MPO from the Appendix in arXiv:2302.14181.
ProductMPO
OptimalFirstOrderEvolutionMPO(BasicTriangularMPO& HamMPO, std::complex<double> Tau)
{
   std::deque<OperatorComponent> Result;

   for (auto const& I : HamMPO)
   {
      // Get the indices for the middle rows and columns.
      std::set<int> MidRows, MidCols;
      for (int i = 1; i < I.Basis1().size()-1; ++i)
         MidRows.insert(i);
      for (int j = 1; j < I.Basis2().size()-1; ++j)
         MidCols.insert(j);

      // Get the Identity, A, B, C and D blocks:
      // (I C D)
      // (0 A B)
      // (0 0 I)
      OperatorComponent Ident = project_columns(project_rows(I, {0}), {0});
      OperatorComponent A = project_columns(project_rows(I, MidRows), MidCols);
      OperatorComponent B = project_columns(project_rows(I, MidRows), {(int) I.Basis2().size()-1});
      OperatorComponent C = project_columns(project_rows(I, {0}), MidCols);
      OperatorComponent D = project_columns(project_rows(I, {0}), {(int) I.Basis2().size()-1});

      // Calculate the four blocks of the optimal first-order evolution MPO.
      OperatorComponent WTopLeft = Ident + Tau * D + 0.5*Tau*Tau * aux_tensor_prod(D, D);
      OperatorComponent WBotLeft = Tau * B + 0.5*Tau*Tau * (aux_tensor_prod(B, D) + aux_tensor_prod(D, B));
      OperatorComponent WTopRight = C + 0.5*Tau * (aux_tensor_prod(C, D) + aux_tensor_prod(D, C));
      OperatorComponent WBotRight = A + 0.5*Tau * (aux_tensor_prod(A, D) + aux_tensor_prod(D, A)
                                                 + aux_tensor_prod(B, C) + aux_tensor_prod(C, B));

      // Construct the evolution MPO.
      OperatorComponent W = tensor_row_sum(tensor_col_sum(WTopLeft, WBotLeft), tensor_col_sum(WTopRight, WBotRight));

      Result.push_back(W);
   }

   return ProductMPO(GenericMPO(Result.begin(), Result.end()));
}

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
      int Verbose = 0;
      int OutputDigits = 0;
      bool Normalize = false;
      bool NoNormalize = false;

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
         ("normalize", prog_opt::bool_switch(&Normalize), "Normalize the wavefunction [default true if timestep is real]")
         ("nonormalize", prog_opt::bool_switch(&NoNormalize), "Don't normalize the wavefunction")
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
         if (Psi.size() < Size)
            std::cerr << "Warning: extending wavefunction unit cell size to " << Size << " sites." << std::endl;
         Psi = repeat(Psi, Size / Psi.size());
         HamMPO = repeat(HamMPO, Size / HamMPO.size());
      }

      //ProductMPO EvolutionMPO = FirstOrderEvolutionMPO(HamMPO, Timestep);
      ProductMPO EvolutionMPO = OptimalFirstOrderEvolutionMPO(HamMPO, Timestep);

      LinearWavefunction PsiL = get_left_canonical(Psi).first;

      if (SaveEvery == 0)
         SaveEvery = N;

      std::cout << "Timestep=" << 0
                << " Time=" << formatting::format_digits(InitialTime, OutputDigits)
                << std::endl;

      double LogAmplitude = Psi.log_amplitude();
      
      for (int tstep = 1; tstep <= N; ++tstep)
      {
         // Apply the evolution MPO.
         auto MI = EvolutionMPO.begin();
         for (auto I = PsiL.begin(); I != PsiL.end(); ++I, ++MI)
            *I = aux_tensor_prod(*MI, *I);

         std::cout << "Timestep=" << tstep
                   << " Time=" << formatting::format_digits(InitialTime+double(tstep)*Timestep, OutputDigits)
                   << std::endl;

         InfiniteWavefunctionLeft PsiSave
            = InfiniteWavefunctionLeft::ConstructPreserveAmplitude(PsiL, Psi.qshift(), Normalize ? 0.0 : LogAmplitude, Verbose);

         //PsiL = get_left_canonical(PsiSave).first;
         LogAmplitude = PsiSave.log_amplitude();

         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
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
