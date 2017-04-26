// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-dmrg.cpp
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


#include "mpo/basic_triangular_mpo.h"
#include "wavefunction/finitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include "common/prog_options.h"
#include <iostream>
#include "common/environment.h"
#include "common/statistics.h"
#include "common/prog_opt_accum.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "mp-algorithms/stateslist.h"
#include "wavefunction/operator_actions.h"

#include "interface/inittemp.h"
#include "mp-algorithms/random_wavefunc.h"

#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/dmrg.h"

namespace prog_opt = boost::program_options;

void SweepRight(DMRG& dmrg, StatesInfo const& SInfo, int SweepNum)
{
   double SweepTruncation = 0;
   dmrg.StartSweep();
   while (dmrg.Site < dmrg.RightStop)
   {
      dmrg.Solve();
      TruncationInfo States = dmrg.TruncateAndShiftRight(SInfo);
      std::cout << "Sweep=" << SweepNum
		<< " Site=" << dmrg.Site
		<< " Energy=";
      if (dmrg.Solver().is_complex())
	 std::cout << dmrg.Solver().LastEnergy();
      else
	 std::cout << dmrg.Solver().LastEnergyReal();
      std::cout << " States=" << States.KeptStates()
		<< " Truncrror=" << States.TruncationError() 
		<< " Iter=" << dmrg.Solver().LastIter()
		<< " Tol=" << dmrg.Solver().LastTol()
		<< '\n';
      SweepTruncation += States.TruncationError();
   }
   std::cout << "Cumumative truncation error for sweep: " << SweepTruncation << '\n';
}

void SweepLeft(DMRG& dmrg, StatesInfo const& SInfo, int SweepNum)
{
   double SweepTruncation = 0;
   dmrg.StartSweep();
   while (dmrg.Site > dmrg.LeftStop)
   {
      dmrg.Solve();
      TruncationInfo States = dmrg.TruncateAndShiftLeft(SInfo);
      std::cout << "Sweep=" << SweepNum
		<< " Site=" << dmrg.Site
		<< " Energy=";
      if (dmrg.Solver().is_complex())
	 std::cout << dmrg.Solver().LastEnergy();
      else
	 std::cout << dmrg.Solver().LastEnergyReal();
      std::cout << " States=" << States.KeptStates()
		<< " Truncrror=" << States.TruncationError()
		<< " Iter=" << dmrg.Solver().LastIter()
		<< " Tol=" << dmrg.Solver().LastTol()
		<< '\n';
      SweepTruncation += States.TruncationError();
   }
   std::cout << "Cumumative truncation error for sweep: " << SweepTruncation << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      int NumIter = 20;
      std::string FName;
      std::string HamStr;
      double FidelityScale = 1.0;
      int MinIter = 4;
      int MaxStates = 100000;
      double MixFactor = 0.01;
      double RandomMixFactor = 0.0;
      bool TwoSite = false;
      int NumSweeps = 10;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      bool NoVariance = false;
      bool UseDGKS = false;
      std::string Solver = "lanczos";
      bool Quiet = false;
      int Verbose = 0;
      double ShiftInvertEnergy = 0.0;
      double MaxTol = 4E-4;  // never use an eigensolver tolerance larger than this
      double MinTol = 1E-16; // lower bound for the eigensolver tolerance - seems we dont really need it
      std::string States = "100";
      double EvolveDelta = 0.0;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
         ("two-site,2", "modify 2 neighboring sites at once (traditional DMRG)")
         ("states,m", prog_opt::value(&States),
          FormatDefault("number of states, or a StatesList", States).c_str())
         ("max-states", prog_opt::value<int>(&MaxStates),
          FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("mix-factor", prog_opt::value(&MixFactor),
          FormatDefault("Mixing coefficient for the density matrix", MixFactor).c_str())
         ("random-mix-factor", prog_opt::value(&RandomMixFactor),
          FormatDefault("Random mixing for the density matrix", RandomMixFactor).c_str())
         ("maxiter", prog_opt::value<int>(&NumIter),
          FormatDefault("Maximum number of Lanczos iterations per step (Krylov subspace size)", NumIter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations per step", MinIter).c_str())
         ("maxtol", prog_opt::value(&MaxTol),
          FormatDefault("Maximum tolerance of the eigensolver", MaxTol).c_str())
         ("sweeps,s", prog_opt::value(&NumSweeps),
          FormatDefault("Number of half-sweeps to perform", NumSweeps).c_str())
         ("Solver,S", prog_opt::value(&Solver),
          FormatDefault("Eigensoler to use ('lanczos', 'arnoldi', 'davidson')", Solver).c_str())
         ("orthogonal", prog_opt::value<std::vector<std::string> >(),
          "force the wavefunction to be orthogonal to this state")
         ("no-variance", prog_opt::bool_switch(&NoVariance), "Don't calculate the variance")
         ("dgks", prog_opt::bool_switch(&UseDGKS), "Use DGKS correction for the orthogonality vectors")
	 ("shift-invert-energy", prog_opt::value(&ShiftInvertEnergy),
	  "For the shift-invert solver, the target energy")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity (can be used more than once)")
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
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      std::cout << "Starting DMRG...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      // Open the wavefunction
      MPWavefunction Wavefunction;
      {
	 pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(FName.c_str(), mp_pheap::CacheSize());
	 Wavefunction = *PsiPtr;
      }
      FiniteWavefunctionLeft Psi = Wavefunction.get<FiniteWavefunctionLeft>();

      // Hamiltonian
      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO;
      
      // get the Hamiltonian from the attributes, if it wasn't supplied
      if (HamStr.empty())
      {
         if (Wavefunction.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: no Hamiltonian specified, use -H option or set wavefunction attribute Hamiltonian.\n";
            return 1;
         }
         HamStr = Wavefunction.Attributes()["Hamiltonian"].as<std::string>();
      }
      else
         Wavefunction.Attributes()["Hamiltonian"] = HamStr;

      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);
      if (HamMPO.size() < Psi.size())
	 HamMPO = repeat(HamMPO, Psi.size() / HamMPO.size());

      // Now we can construct the actual DMRG object
      DMRG dmrg(Psi, HamMPO, Verbose);

      dmrg.UseDGKS = UseDGKS;
      dmrg.Solver().SetSolver(Solver);

      dmrg.Solver().MaxTol = MaxTol;
      dmrg.Solver().MinTol = MinTol;
      dmrg.Solver().MinIter = MinIter;
      dmrg.Solver().MaxIter = NumIter;
      dmrg.Solver().FidelityScale = FidelityScale;
      dmrg.Solver().Verbose = Verbose;
      dmrg.Solver().EvolveDelta = EvolveDelta;
      dmrg.Solver().SetShiftInvertEnergy(ShiftInvertEnergy);

      dmrg.MixingInfo.MixFactor = MixFactor;
      dmrg.MixingInfo.RandomMixFactor = RandomMixFactor;

      StatesInfo SInfo;
      SInfo.MinStates = 1;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      std::cout << SInfo << '\n';

      StatesList MyStates(States);
      if (vm.count("steps") && MyStates.size() == 1)
      {
         MyStates.Repeat(NumSweeps);
      }
      std::cout << MyStates << '\n';

      std::cout << "Density matrix mixing coefficient: " << MixFactor << std::endl;
      std::cout << "Number of Lanczos iterations: " << NumIter << std::endl;
      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      std::cout << "Using solver: " << Solver << std::endl;

      for (int Sweeps = 0; Sweeps < MyStates.size(); ++Sweeps)
      {
	 SInfo.MaxStates = MyStates[Sweeps].NumStates;
         if (Sweeps % 2 == 0)
            SweepLeft(dmrg, SInfo, Sweeps+1);
         else
            SweepRight(dmrg, SInfo, Sweeps+1);

#if 0
         // the dmrg.Wavefunction() is not normalized anymore
         double Norm2 = norm_2_sq(dmrg.Wavefunction());

         // We need to re-calculate the energy, since it will have changed slightly after the truncation
         double E = dmrg.Energy()/Norm2;
         std::cout << "E = " << E << '\n';

         if (CalculateH2)
         {
            double h2 = std::abs(expectation(dmrg.Wavefunction(), Ham2, dmrg.Wavefunction()))/Norm2;
            double nh2 = h2 - E*E;
            std::cout << "(H-E)^2 = " << nh2 << '\n';
         }

         Psi = dmrg.Wavefunction();
         double Overlap = dmrg.FidelityLoss();

         std::cout << "Wavefunction difference from last half-sweep = " << Overlap << '\n';
         OldPsi = Psi;
#endif
      }

      // finished the iterations.
      std::cout << "Orthogonalizing wavefunction...\n";
      Wavefunction.Wavefunction() = dmrg.Wavefunction();

      // any other attributes?
      Wavefunction.Attributes()["LastEnergy"] = dmrg.Solver().LastEnergy();
      Wavefunction.SetDefaultAttributes();

      // History log
      Wavefunction.AppendHistory(EscapeCommandline(argc, argv));

      // save wavefunction
      pvalue_ptr<MPWavefunction> P(new MPWavefunction(Wavefunction));
      pheap::ShutdownPersistent(P);

      ProcControl::Shutdown();
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
