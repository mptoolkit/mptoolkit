// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-itebd.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mp-algorithms/tebd.h"
#include "mpo/basic_triangular_mpo.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "parser/number-parser.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/prog_options.h"
#include "common/statistics.h"
#include "lattice/infinite-parser.h"
#include "common/formatting.h"
#include <cctype>
#include "interface/inittemp.h"
#include "common/openmp.h"

namespace prog_opt = boost::program_options;

void DoEvenSlice(std::deque<StateComponent>& Psi,
                 std::deque<RealDiagonalOperator>& Lambda,
                 std::vector<SimpleOperator> const& UEven,
                 StatesInfo const& SInfo,
                 int Verbose)
{
   // physical state is A B (Lambda) C D (Lambda) ...
   // All A,B,C,D,... are left canonical
   // After even slice, the state is
   // A (Lambda) B C (Lambda) D E (Lambda) ...
   unsigned Sz = Psi.size();
   int MaxStates = 0;
   #pragma omp parallel for shared(Psi, Lambda, UEven, SInfo)
   for (unsigned i = 0; i < Sz; i += 2)
   {
      TruncationInfo Info = DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UEven[i/2], SInfo);
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
   std::cout << "Even slice finished, max states=" << MaxStates << '\n';
   std::cout << std::flush;
}

void DoOddSlice(std::deque<StateComponent>& Psi,
                std::deque<RealDiagonalOperator>& Lambda,
                QuantumNumber const& QShift,
                std::vector<SimpleOperator> const& UOdd,
                StatesInfo const& SInfo,
                int Verbose)
{
   // In preparation for the odd slice, we rotate to
   // Z A (Lambda) B C (Lambda) D E (Lambda) ...
   // After the odd slice, we have
   // Z (Lamda) A B (Lambda) C D (Lambda) ...
   // which we need to rotate back to
   // A B (Lambda) C D (Lambda) ..... Z (Lambda)
   unsigned Sz = Psi.size();
   Psi.push_front(delta_shift(Psi.back(), QShift));
   Psi.pop_back();
   int MaxStates = 0;
   #pragma omp parallel for shared(Psi, Lambda, UOdd, SInfo)
   for (unsigned i = 0; i < Sz; i += 2)
   {
      TruncationInfo Info = DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UOdd[i/2], SInfo);
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
   std::cout << "Odd slice finished, max states=" << MaxStates << '\n';
   std::cout << std::flush;

   Psi.push_back(delta_shift(Psi.front(), adjoint(QShift)));
   Psi.pop_front();
   Lambda.push_back(delta_shift(Lambda.front(), adjoint(QShift)));
   Lambda.pop_front();
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
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      int OutputDigits = 0;
      int Coarsegrain = 1;
      std::string DecompositionStr = "optimized4-11";

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "operator to use for the Hamiltonian (wavefunction attribute \"EvolutionHamiltonian\") [required]")
	 ("wavefunction,w", prog_opt::value(&InputFile), "input wavefunction [required]")
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
         ("coarsegrain", prog_opt::value(&Coarsegrain),
          "coarse-grain N-to-1 sites")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1 || (vm.count("timestep") < 1 && vm.count("betastep") < 1))
      {
         print_copyright(std::cerr, "tools", "mp-itebd");
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "\nAvailable decompositions:\n";
         for (auto const& d : Decompositions)
         {
            std::cerr << d.first << " : " << d.second.Description_ << '\n';
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
      if (decomp.Order_ == 0)
      {
         std::cerr << "mp-itebd: fatal: invalid decomposition\n";
         return 1;
      }

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

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

      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

      //Coarse-grain, if necessary
      if (vm.count("coarsegrain"))
      {
         HamMPO = coarse_grain(HamMPO, Coarsegrain);
      }

      // Make sure the unit cell is an even size
      int UnitCellSize = statistics::lcm(Psi.size(), HamMPO.size(), 2);

      // and adjust the wavefunction and hamiltonian to match the unit cell size
      if (Psi.size() != UnitCellSize)
      {
         std::cerr << "mp-itebd: warning: extending wavefunction unit cell to " << UnitCellSize << " sites.\n";
         Psi = repeat(Psi, UnitCellSize / Psi.size());
      }
      HamMPO = repeat(HamMPO, UnitCellSize / HamMPO.size());

      // Assemble the Hamiltonian into the bond terms
      std::vector<SimpleOperator> BondH(UnitCellSize);
      auto SplitTerms = SplitMPO(HamMPO);
      for (int i = 0; i < SplitTerms.size(); ++i)
      {
         for (auto const& op : SplitTerms[i])
         {
            if (op.size() == 1)
            {
               // split single-site operators over the two bonds, if i=0 then
               // the split wraps around the unit cell
               SimpleRedOperator p = tensor_prod(op[0](0,0),HamMPO[(i+1)%UnitCellSize](0,0));
               BondH[i] += 0.5 * p.scalar();
               SimpleRedOperator q = tensor_prod(HamMPO[(i+UnitCellSize-1)%UnitCellSize](0,0), op[0](0,0));
               BondH[(i+UnitCellSize-1)%UnitCellSize] += 0.5 * q.scalar();
            }
            else if (op.size() == 2)
            {
                  BondH[i] += coarse_grain(op).scalar();
            }
            else
            {
               throw std::runtime_error("fatal: operator has support over more than 2 sites.");
            }
         }
      }

      // Exponentiate the bond operators and split into even and odd slices
      std::vector<std::vector<SimpleOperator>> EvenU;
      for (auto x : decomp.a_)
      {
         std::vector<SimpleOperator> Terms;
         for (int i = 0; i < BondH.size(); i += 2)
         {
            Terms.push_back(Exponentiate(-Timestep*std::complex<double>(0,x) * BondH[i]));
         }
         EvenU.push_back(std::move(Terms));
      }

      std::vector<std::vector<SimpleOperator>> OddU;
      // The odd slice uses a rotation of the unit cell that takes the right-most site and puts it on the left.
      // This means that the first odd-bond operator we need is the one that wraps around, which is at the end of the list.
      // To allow for this, we reserve the first element of the Terms vector and move the last element to the front
      for (auto x : decomp.b_)
      {
         std::vector<SimpleOperator> Terms;
         Terms.push_back(Exponentiate(-Timestep*std::complex<double>(0,x) * BondH[BondH.size()-2]));
         for (int i = 1; i < BondH.size()-2; i += 2)
         {
            Terms.push_back(Exponentiate(-Timestep*std::complex<double>(0,x) * BondH[i]));
         }
         OddU.push_back(std::move(Terms));
      }

      // If we have an odd number of terms (the usual case), then we can wrap around the last even slice
      // if we are continuing the evolution beyond the current timestep.  This is the sum of a[last] + a[0] terms
      std::vector<SimpleOperator> EvenContinuation;
      if (decomp.a_.size() == decomp.b_.size()+1)
      {
         double x = decomp.a_.front() + decomp.a_.back();
         for (int i = 0; i < BondH.size(); i += 2)
         {
            EvenContinuation.push_back(Exponentiate(-Timestep*std::complex<double>(0,x) * BondH[i]));
         }
      }

      std::cout << "Using decomposition " << DecompositionStr << '\n';
      if (Verbose > 0)
      {
         std::cout << "Number of even slices: " << EvenU.size() << '\n';
         std::cout << "Number of odd slices: " << OddU.size() << '\n';
      }

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = States;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      std::cout << SInfo << '\n';

      QuantumNumber QShift = Psi.qshift();

      std::deque<StateComponent> PsiVec(Psi.begin(), Psi.end());
      std::deque<RealDiagonalOperator> Lambda;

      for (int i = 2; i <= Psi.size(); i += 2)
      {
         Lambda.push_back(Psi.lambda(i));
      }

      if (SaveEvery == 0)
         SaveEvery = N;

      // If Continue then we merge the final (even) slice of one timestep with the first slice of
      // the next timestep
      bool Continue = false;
      int tstep = 0;

      while (tstep < N)
      {
         if (Continue)
         {
            if (Verbose > 1)
            {
               std::cout << "Merge slice\n";
            }
            DoEvenSlice(PsiVec, Lambda, EvenContinuation, SInfo, Verbose);
         }
         else
         {
            DoEvenSlice(PsiVec, Lambda, EvenU[0], SInfo, Verbose);
         }

         DoOddSlice(PsiVec, Lambda, QShift, OddU[0], SInfo, Verbose);
         for (int bi = 1; bi < OddU.size(); ++bi)
         {
            DoEvenSlice(PsiVec, Lambda, EvenU[bi], SInfo, Verbose);
            DoOddSlice(PsiVec, Lambda, QShift, OddU[bi], SInfo, Verbose);
         }

         // do we do a continuation?
         Continue = EvenU.size() > OddU.size();

         ++tstep;
         std::cout << "Timestep " << formatting::format_complex(tstep)
                   << " time " << formatting::format_complex(InitialTime+double(tstep)*Timestep) << '\n';

         // do we save the wavefunction?
         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            LinearWavefunction Psi;
            if (EvenU.size() > OddU.size())
            {
               CHECK_EQUAL(EvenU.size(), OddU.size()+1);
               std::cout << "Doing final slice before saving wavefunction.\n";
               // do the final slice to finish the timstep.  Make a copy of the wavefunction since it is better to
               // avoid a truncation step and 'continue' the wavefunction by wrapping around the next timestep.
               std::deque<StateComponent> PsiVecSave = PsiVec;
               std::deque<RealDiagonalOperator> LambdaSave = Lambda;
               DoEvenSlice(PsiVecSave, LambdaSave, EvenU.back(), SInfo, Verbose);
               Psi = LinearWavefunction::FromContainer(PsiVecSave.begin(), PsiVecSave.end());
            }
            else
               Psi = LinearWavefunction::FromContainer(PsiVec.begin(), PsiVec.end());

            // save the wavefunction
            std::cout << "Saving wavefunction\n";
            MPWavefunction Wavefunction;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), OutputDigits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), OutputDigits);
            InfiniteWavefunctionLeft PsiL = InfiniteWavefunctionLeft::Construct(Psi, QShift);
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
