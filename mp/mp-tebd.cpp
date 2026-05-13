// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-tebd.cpp
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp-algorithms/tebd.h"
#include "mp-algorithms/time-dependent-mpo.h"
#include "mpo/basic_triangular_mpo.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "parser/number-parser.h"
#include <iostream>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/prog_options.h"
#include "lattice/infinite-parser.h"
#include "interface/inittemp.h"
#include "common/openmp.h"
#include "common/formatting.h"
#include "tensor/reducible.h"
#include "parser/parser.h"
#include <cctype>

namespace prog_opt = boost::program_options;

void DoEvenSlice(std::deque<StateComponent>& Psi,
                 std::deque<RealDiagonalOperator>& Lambda,
                 double& LogAmplitude,
                 std::vector<SimpleOperator> const& UEven,
                 StatesInfo const& SInfo,
                 int Verbose)
{
   // physical state is A B (Lambda) C D (Lambda) ...
   // All A,B,C,D,... are left canonical
   // After even slice, the state is
   // A (Lambda) B C (Lambda) D E (Lambda) ...
   // If the number of sites is even, we replicate the last Lambda matrix, so that
   // we always have a final Lambda matrix at the end
   unsigned Sz = Psi.size();
   if (Sz%2 == 0)
   {
      Lambda.push_back(Lambda.back());
   }
   int MaxStates = 0;
   double DeltaLogAmplitude = 0;
   #pragma omp parallel for reduction(+:DeltaLogAmplitude)
   for (unsigned i = 0; i < Sz-1; i += 2)
   {
      TruncationInfo Info = DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], DeltaLogAmplitude, UEven[i/2], SInfo);
      if (Verbose > 0)
      {
         std::cout << "Bond=" << (i+1)
                   << " Entropy=" << Info.TotalEntropy()
                   << " States=" << Info.KeptStates()
                   << " Trunc=" << Info.TruncationError()
                   << '\n';
      }
      #pragma omp critical
      MaxStates = std::max(MaxStates, Info.KeptStates());
   }
   LogAmplitude += DeltaLogAmplitude;
   std::cout << "Even slice finished, max states=" << MaxStates << '\n';
   std::cout << std::flush;
}

void DoOddSlice(std::deque<StateComponent>& Psi,
                std::deque<RealDiagonalOperator>& Lambda,
                double& LogAmplitude,
                std::vector<SimpleOperator> const& UOdd,
                StatesInfo const& SInfo,
                int Verbose)
{
   // Physical state is
   // A (Lambda) B C (Lambda) D E (Lambda) ... (all left-canonical)
   // After an odd slice, the physical state is
   // A B (Lambda) C D (Lambda) ...
   // Each TEBD iteration turns a set B C (Lambda) -> B (Lambda) C
   // The first (Lambda) matrix is unused; we remove it from the container.
   // If the number of sites is odd, we copy the last Lambda matrix so that we always have a
   // 1x1 identity lambda at the end

   // If the system size is odd:
   // Physical state is
   // A (Lambda) B C (Lambda)
   // We want to access lambda matrices at indices 1,3,5,...
   // after we remove the first entry, we index lambda 0,2,4,...
   // on exit the first lambda in the container comes from psi[1] psi[2]
   // After odd slice, physical state is
   // A B (Lambda) C

   unsigned Sz = Psi.size();
   if (Sz%2 == 1)
   {
      Lambda.push_back(Lambda.back());
   }
   Lambda.pop_front();

   int MaxStates = 0;
   double DeltaLogAmplitude = 0;
   #pragma omp parallel for reduction(+:DeltaLogAmplitude)
   for (unsigned i = 1; i < Sz-1; i += 2)
   {
      TruncationInfo Info = DoTEBD(Psi[i], Psi[i+1], Lambda[(i-1)/2], DeltaLogAmplitude, UOdd[(i-1)/2], SInfo);
      if (Verbose > 0)
      {
         std::cout << "Bond=" << i
                   << " Entropy=" << Info.TotalEntropy()
                   << " States=" << Info.KeptStates()
                   << " Trunc=" << Info.TruncationError()
                   << '\n';
      }
      #pragma omp critical
         MaxStates = std::max(MaxStates, Info.KeptStates());
   }
   LogAmplitude += DeltaLogAmplitude;
   std::cout << "Odd slice finished, max states=" << MaxStates << '\n';
   std::cout << std::flush;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string TimestepStr;
      std::string BetastepStr;
      std::string InitialTimeStr;
      std::string InitialBetaStr;
      int N = 1;
      int SaveEvery = 1;
      std::string OpStr;
      std::string InputFile;
      std::string OutputPrefix;
      std::string HamStr;
      int MinStates = 1;
      int States = 100000;
      bool Normalize = false;
      bool NoNormalize = false;
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      int OutputDigits = 0;
      int Coarsegrain = 1;
      int MagnusOrder = 2;
      int MagnusQuadrature = 0;
      std::string TimeVar = "t";
      std::string DecompositionStr = "optimized4-11";
      bool TimeDependent = false;
      bool Quiet = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "operator to use for the Hamiltonian (wavefunction attribute \"EvolutionHamiltonian\")")
         ("wavefunction,w", prog_opt::value(&InputFile), "input wavefunction")
         ("output,o", prog_opt::value(&OutputPrefix), "prefix for saving output files")
         ("timestep,t", prog_opt::value(&TimestepStr), "timestep (required)")
         ("betastep,b", prog_opt::value(&BetastepStr), "betastep (alternative to timestep)")
         ("decomposition,c", prog_opt::value(&DecompositionStr), FormatDefault("choice of decomposition", DecompositionStr).c_str())
         ("num-timesteps,n", prog_opt::value(&N), FormatDefault("number of timesteps to calculate", N).c_str())
         ("save-timesteps,s", prog_opt::value(&SaveEvery), "save the wavefunction every s timesteps")
         ("min-states", prog_opt::value<int>(&MinStates),
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("states,m", prog_opt::value(&States),
          FormatDefault("number of states", States).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("magnus", prog_opt::value(&MagnusOrder), FormatDefault("For time-dependent Hamiltonians, use this order of the Magnus expansion", MagnusOrder).c_str())
         ("magnus-quadrature", prog_opt::value(&MagnusQuadrature),
          FormatDefault("Gauss-Legendre quadrature order for time-dependent Magnus terms (0 selects the default for --magnus)", MagnusQuadrature).c_str())
         ("timevar", prog_opt::value(&TimeVar), FormatDefault("The time variable for time-dependent Hamiltonians", TimeVar).c_str())
         ("coarsegrain", prog_opt::value(&Coarsegrain),
          "coarse-grain N-to-1 sites")
         ("normalize", prog_opt::bool_switch(&Normalize), "Normalize the wavefunction [default if timestep is real]")
         ("nonormalize", prog_opt::bool_switch(&NoNormalize), "Don't normalize the wavefunction [default if the timestep is complex]")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1
          ||  (vm.count("timestep") < 1 && vm.count("betastep") < 1))
      {
         print_copyright(std::cerr, "tools", "mp-tebd");
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Available decompositions:\n";
         for (auto const& d : Decompositions)
         {
            std::cerr << d.first << " : " << d.second.description() << '\n';
         }
         return 1;
      }

      omp::initialize(Verbose);
      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      LTSDecomposition decomp;
      for (auto const& d : Decompositions)
      {
         if (d.first == DecompositionStr)
            decomp = d.second;
      }
      if (decomp.order() == 0)
      {
         std::cerr << "mp-tebd: fatal: invalid decomposition\n";
         return 1;
      }

      try
      {
         ResolveMagnusQuadratureOrder(MagnusOrder, MagnusQuadrature);
      }
      catch (std::exception const& e)
      {
         std::cerr << "mp-tebd: fatal: " << e.what() << '\n';
         return 1;
      }

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);

      FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();

      // Ooutput prefix - if it wasn't specified then use the wavefunction attribute, or
      // fallback to the wavefunction name
      if (OutputPrefix.empty())
         OutputPrefix = PsiPtr->Attributes()["Prefix"].as<std::string>();

      if (OutputPrefix.empty())
         OutputPrefix = InputFile;

      // Get the initial time & beta from the wavefunction attributes
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
         Timestep += std::complex<double>(0.0,-1.0)* ParseNumber(BetastepStr);

      if (Timestep.imag() == 0.0)
         Normalize = true;
      if (NoNormalize)
         Normalize = false;

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      std::cout << "Starting TEBD.\nTrotter-Suziki decomposition is " << DecompositionStr << "(" << decomp.description() << ")\n"
                << "Timestep = " << formatting::format_complex(Timestep) << '\n';

      if (OutputDigits == 0)
      {
         OutputDigits = std::max(formatting::digits(Timestep), formatting::digits(InitialTime));
      }

      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;

      // get the Hamiltonian from the attributes, if it wasn't supplied.
      // Get it from the EvolutionHamiltonian, if it exists,
      // or from Hamiltonian
      if (HamStr.empty())
      {
         if (PsiPtr->Attributes().count("EvolutionHamiltonian") == 0)
         {
            if (PsiPtr->Attributes().count("Hamiltonian") == 0)
            {
               std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute EvolutionHamiltonian.\n";
               return 1;
            }
            HamStr = PsiPtr->Attributes()["Hamiltonian"].as<std::string>();
         }
         else
         {
            HamStr = PsiPtr->Attributes()["EvolutionHamiltonian"].as<std::string>();
         }
      }

      std::string HamOperator;
      std::tie(HamOperator, Lattice) = ParseOperatorStringAndLattice(HamStr);

      // Attempt to convert the operator into an MPO.  If it works, then the
      // Hamiltonian is time-independent.  If it fails, assume it failed because
      // the time variable is not defined yet; other operator errors will be
      // reported when the time-dependent MPO is built.
      try
      {
         HamMPO = ParseTriangularOperator(Lattice, HamOperator);
         HamMPO = PrepareFiniteTEBDHamiltonian(HamMPO, Psi.size());
      }
      catch (Parser::ParserError& e)
      {
         if (Verbose > 1)
         {
            std::cerr << "Parser error converting the Hamiltonian to an MPO - assuming the Hamiltonian is time-dependent.\n";
         }
         TimeDependent = true;
      }
      catch (...)
      {
         throw;
      }

      TEBDHamiltonianGates Gates;
      if (!TimeDependent)
      {
         Gates = AssembleFiniteTEBDHamiltonian(HamMPO, Timestep, decomp);
      }

      std::cout << "Using decomposition " << DecompositionStr << '\n';
      if (Verbose > 0)
      {
         if (!TimeDependent)
         {
            std::cout << "Number of even slices: " << Gates.EvenU.size() << '\n';
            std::cout << "Number of odd slices: " << Gates.OddU.size() << '\n';
         }
      }

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = States;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      std::cout << SInfo << '\n';

      double LogAmplitude = 0.0; // normalization is contained in the lambda matrices

      std::deque<StateComponent> PsiVec(Psi.begin(), Psi.end());
      std::deque<RealDiagonalOperator> Lambda;

      // The Lambda matrices we need are on the right hand side of every pair of A-matrices, so sites
      // 2,4,6,...
      // But we always want the edge lambda, so if size is odd, we add one additional lambda matrix. This
      // is to cover the case where the system size is odd, and when we do a odd slice, we need this matrix.
      for (int i = 2; i < Psi.size(); i += 2)
      {
         Lambda.push_back(Psi.lambda(i));
      }
      Lambda.push_back(Psi.lambda_r());

      if (SaveEvery == 0)
         SaveEvery = N;

      // If Continue then we merge the final (even) slice of one timestep with the first slice of
      // the next timestep
      bool Continue = false;
      int tstep = 0;

      while (tstep < N)
      {
         // If the Hamiltonian is time-dependent, build the gates over the
         // actual time interval covered by each Suzuki-Trotter slice.
         if (TimeDependent)
         {
            Gates = AssembleFiniteTimeDependentTEBDHamiltonian(Lattice, HamOperator, TimeVar,
                                                               InitialTime + double(tstep)*Timestep, Timestep,
                                                               decomp, MagnusOrder, MagnusQuadrature,
                                                               Psi.size());
         }

         if (Continue)
         {
            if (Verbose > 1)
            {
               std::cout << "Merge slice\n";
            }
            DoEvenSlice(PsiVec, Lambda, LogAmplitude, Gates.EvenContinuation, SInfo, Verbose);
         }
         else
         {
            DoEvenSlice(PsiVec, Lambda, LogAmplitude, Gates.EvenU[0], SInfo, Verbose);
         }

         DoOddSlice(PsiVec, Lambda, LogAmplitude, Gates.OddU[0], SInfo, Verbose);
         for (int bi = 1; bi < Gates.OddU.size(); ++bi)
         {
            DoEvenSlice(PsiVec, Lambda, LogAmplitude, Gates.EvenU[bi], SInfo, Verbose);
            DoOddSlice(PsiVec, Lambda, LogAmplitude, Gates.OddU[bi], SInfo, Verbose);
         }

         // do we do a continuation?
         Continue = Gates.EvenU.size() > Gates.OddU.size();

         ++tstep;
         std::cout << "Timestep " << formatting::format_complex(tstep)
                   << " time " << formatting::format_complex(InitialTime+double(tstep)*Timestep) << '\n';

         // do we save the wavefunction?
         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            LinearWavefunction Psi;
            double ThisLogAmplitude = LogAmplitude;
            if (Gates.EvenU.size() > Gates.OddU.size())
            {
               CHECK_EQUAL(Gates.EvenU.size(), Gates.OddU.size()+1);
               std::cout << "Doing final slice before saving wavefunction.\n";
               // do the final slice to finish the timstep.  Make a copy of the wavefunction since it is better to
               // avoid a truncation step and 'continue' the wavefunction by wrapping around the next timestep.
               std::deque<StateComponent> PsiVecSave = PsiVec;
               std::deque<RealDiagonalOperator> LambdaSave = Lambda;
               DoEvenSlice(PsiVecSave, LambdaSave, ThisLogAmplitude, Gates.EvenU.back(), SInfo, Verbose);
               Psi = LinearWavefunction::FromContainer(PsiVecSave.begin(), PsiVecSave.end());
            }
            else
               Psi = LinearWavefunction::FromContainer(PsiVec.begin(), PsiVec.end());

            // save the wavefunction
            std::cout << "Saving wavefunction\n";
            MPWavefunction Wavefunction;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), OutputDigits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), OutputDigits);
            FiniteWavefunctionLeft PsiL = FiniteWavefunctionLeft::Construct(Psi);
            if (!Normalize)
            {
               PsiL *= std::exp(ThisLogAmplitude);
            }
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
