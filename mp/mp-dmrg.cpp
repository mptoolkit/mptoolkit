// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-dmrg.cpp
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
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
#include "common/formatting.h"
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

bool Bench = (getenv_or_default("MP_BENCHFILE", "") != std::string());
std::ofstream BenchFile(getenv_or_default("MP_BENCHFILE", ""), std::ios_base::out | std::ios_base::trunc);

bool Flush = false;  // set to true to flush standard output every step

struct ExpansionInfo
{
   double IncrementFactor;        // If the number of states is increasing from m to m', then add IncrementFactor*(m'-m) states
   double ExpandFactor;           // Add additional ExpandFactor * m' states
   int ExpandPerSector;           // Add this many states per available quantum number sector

   bool should_expand() const
   {
      return IncrementFactor > 0.0 || ExpandFactor > 0.0 || ExpandPerSector > 0;
   }

   ExpansionInfo()
      : IncrementFactor(0.0), ExpandFactor(0.0), ExpandPerSector(0) {}
};

ExpansionInfo operator*(ExpansionInfo e, double factor)
{
   e.ExpandFactor *= factor;
   e.ExpandPerSector = int(std::round(e.ExpandPerSector * factor));
   return e;
}

std::ostream& operator<<(std::ostream& out, ExpansionInfo const& Info)
{
   out << "increment: " << Info.IncrementFactor << " expansion factor: " << Info.ExpandFactor << " per sector: " << Info.ExpandPerSector;
   return out;
}

void SweepRight(DMRG& dmrg, int SweepNum, StatesInfo const& States, ExpansionInfo const& PreExpand, ExpansionInfo const& PostExpand, int NumStatesKeepNext)
{
   double SweepTruncation = 0;
   dmrg.StartSweep();
   while (dmrg.Site < dmrg.RightStop)
   {
      dmrg.StartIteration();
      int CurrentEnvStates = dmrg.BasisTotalDimension2();
      int DesiredStates = States.MaxStates;
      int ExtraStates = int(std::ceil(PreExpand.IncrementFactor*std::max(DesiredStates-CurrentEnvStates, 0) + PreExpand.ExpandFactor*DesiredStates));

      if ((ExtraStates > 0 || PreExpand.ExpandPerSector > 0)) // && dmrg.Site < dmrg.RightStop)
      {
         CurrentEnvStates = dmrg.ExpandRightEnvironment(CurrentEnvStates+ExtraStates, PreExpand.ExpandPerSector);
      }
      dmrg.Solve();
      int Delta = int(std::ceil(PostExpand.IncrementFactor*std::max(NumStatesKeepNext-DesiredStates, 0) + PostExpand.ExpandFactor*States.MaxStates));
      //TRACE(Delta)(PostExpand.IncrementFactor)(PostExpand.ExpandFactor);
      TruncationInfo Info = dmrg.TruncateAndShiftRight(States, Delta, PostExpand.ExpandPerSector);
      std::cout << "Sweep=" << SweepNum
         << " Site=" << dmrg.Site
         << " Energy=" << formatting::format_complex(dmrg.Solver().LastEnergy())
         << " Env=" << CurrentEnvStates
         << " States=" << Info.KeptStates()
         << " Extra=" << Info.ExtraStates()
         << " Truncrror=" << Info.TruncationError()
         << " FidelityLoss=" << dmrg.Solver().LastFidelityLoss()
         << " Iter=" << dmrg.Solver().LastIter()
         << " Tol=" << dmrg.Solver().LastTol()
         << '\n';
      if (Flush)
         std::cout << std::flush;
      SweepTruncation += Info.TruncationError();
      //dmrg.EndIteration();
      if (Bench)
         BenchFile << ProcControl::GetElapsedTime() << ' ' << SweepNum << ' ' << dmrg.Site << ' ' << Info.KeptStates() << ' ' << formatting::format_complex(dmrg.Solver().LastEnergy()) << ' ' << Info.TruncationError() << ' ' << dmrg.Solver().LastFidelityLoss() << ' ' << dmrg.Solver().LastIter() << ' ' << dmrg.Solver().LastTol() << '\n';
   }
   std::cout << "Cumumative truncation error for sweep: " << SweepTruncation << '\n';
}

void SweepLeft(DMRG& dmrg, int SweepNum, StatesInfo const& States, ExpansionInfo const& PreExpand, ExpansionInfo const& PostExpand, int NumStatesKeepNext)
{
   double SweepTruncation = 0;
   dmrg.StartSweep();
   while (dmrg.Site > dmrg.LeftStop)
   {
      dmrg.StartIteration();
      int CurrentEnvStates = dmrg.BasisTotalDimension1();
      int DesiredStates = States.MaxStates;
      int ExtraStates = int(std::ceil(PreExpand.IncrementFactor*std::max(DesiredStates-CurrentEnvStates, 0) + PreExpand.ExpandFactor*DesiredStates));

      if ((ExtraStates > 0 || PreExpand.ExpandPerSector > 0)) // && dmrg.Site < dmrg.RightStop)
      {
         CurrentEnvStates = dmrg.ExpandLeftEnvironment(CurrentEnvStates+ExtraStates, PreExpand.ExpandPerSector);
      }
      dmrg.Solve();
      int Delta = int(std::ceil(PostExpand.IncrementFactor*std::max(NumStatesKeepNext-DesiredStates, 0) + PostExpand.ExpandFactor*States.MaxStates));
      TruncationInfo Info = dmrg.TruncateAndShiftLeft(States, Delta, PostExpand.ExpandPerSector);
      std::cout << "Sweep=" << SweepNum
         << " Site=" << dmrg.Site
         << " Energy=" << formatting::format_complex(dmrg.Solver().LastEnergy())
         << " Env=" << CurrentEnvStates
         << " States=" << Info.KeptStates()
         << " Extra=" << Info.ExtraStates()
         << " Truncrror=" << Info.TruncationError()
         << " FidelityLoss=" << dmrg.Solver().LastFidelityLoss()
         << " Iter=" << dmrg.Solver().LastIter()
         << " Tol=" << dmrg.Solver().LastTol()
         << '\n';
      if (Flush)
         std::cout << std::flush;
      SweepTruncation += Info.TruncationError();
      //dmrg.EndIteration();
      if (Bench)
         BenchFile << ProcControl::GetElapsedTime() << ' ' << SweepNum << ' ' << dmrg.Site << ' ' << Info.KeptStates() << ' ' << formatting::format_complex(dmrg.Solver().LastEnergy()) << ' ' << Info.TruncationError() << ' ' << dmrg.Solver().LastFidelityLoss() << ' ' << dmrg.Solver().LastIter() << ' ' << dmrg.Solver().LastTol() << '\n';
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
      bool TwoSite = false;
      int NumSweeps = 10;
      double TruncCutoff = -1;
      double EigenCutoff = 1e-30;
      int SubspaceSize = 30;
      bool UsePreconditioning = false;
      bool UseDGKS = false;
      std::string Solver = "lanczos";
      bool Quiet = false;
      int Verbose = 0;
      std::complex<double> ShiftInvertEnergy = 0.0;
      double MaxTol = 4E-4;  // never use an eigensolver tolerance larger than this
      double MinTol = 1E-16; // lower bound for the eigensolver tolerance - seems we dont really need it
      std::string StatesString = "100";
      double EvolveDelta = 0.0;
      ExpansionInfo PreExpand;
      ExpansionInfo PostExpand;
      double InitialFidelity = 1e-7;
      std::string PreExpandAlgo = PreExpansionAlgorithm().Name();
      std::string PostExpandAlgo = PostExpansionAlgorithm().Name();
      bool NoGreedy = false;  // set to false to expand the basis quickly, keeping enough states for the folllowing sweep
      bool ProjectTwoSiteTangent = false;
      OversamplingInfo Oversampling(10, 1.0, 1); // Oversample extra state per sector

      // Defaults for expansion
      PreExpand.IncrementFactor = 1.0;
      PreExpand.ExpandFactor = 0.1;
      PreExpand.ExpandPerSector = 0;

      PostExpand.IncrementFactor = 1.0;
      PostExpand.ExpandFactor = 0.05;
      PostExpand.ExpandPerSector = 0;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value(&FName),
          "wavefunction to apply DMRG (required)")
         ("states,m", prog_opt::value(&StatesString),
          FormatDefault("number of states, or a StatesList", StatesString).c_str())
         ("max-states", prog_opt::value<int>(&MaxStates),
          FormatDefault("Maximum number of states to keep", MaxStates).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
          ("pre", prog_opt::value(&PreExpandAlgo), FormatDefault("Pre-expansion algorithm, choices are " + PreExpansionAlgorithm::ListAvailable(), PreExpandAlgo).c_str())
          ("pre-increment", prog_opt::value(&PreExpand.IncrementFactor), FormatDefault("Pre-expansion growth factor for basis size increase", PreExpand.IncrementFactor).c_str())
         ("pre-factor", prog_opt::value(&PreExpand.ExpandFactor), FormatDefault("Pre-expansion factor", PreExpand.ExpandFactor).c_str())
         ("pre-sector", prog_opt::value(&PreExpand.ExpandPerSector), "Pre-expansion number of additional states in each quantum number sector [default 0 for fullsvd, rsvd; default 1 for range, random]")
         ("post", prog_opt::value(&PostExpandAlgo), FormatDefault("Post-expansion algorithm, choices are " + PostExpansionAlgorithm::ListAvailable(), PostExpandAlgo).c_str())
         ("post-increment", prog_opt::value(&PostExpand.IncrementFactor), FormatDefault("Post-expansion growth factor for basis size increase", PostExpand.IncrementFactor).c_str())
         ("post-factor", prog_opt::value(&PostExpand.ExpandFactor), FormatDefault("Post-expansion factor", PostExpand.ExpandFactor).c_str())
         ("post-sector", prog_opt::value(&PostExpand.ExpandPerSector), "Post-expansion number of additional states in each quantum number sector [default 0 for fullsvd, rsvd; default 1 for range, random]")
         ("nogreedy", prog_opt::bool_switch(&NoGreedy), "Don't expand the basis one sweep ahead")
         ("twositetangent", prog_opt::bool_switch(&ProjectTwoSiteTangent), "Project onto the two-site tangent space during pre-expansion")
         ("evolve", prog_opt::value(&EvolveDelta),
          "Instead of Lanczos, do imaginary time evolution with this timestep")
         ("maxiter", prog_opt::value<int>(&NumIter),
          FormatDefault("Maximum number of Lanczos iterations per step (Krylov subspace size)", NumIter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations per step", MinIter).c_str())
         ("maxtol", prog_opt::value(&MaxTol),
          FormatDefault("Maximum tolerance of the eigensolver", MaxTol).c_str())
         ("sweeps,s", prog_opt::value(&NumSweeps),
          FormatDefault("Number of half-sweeps to perform", NumSweeps).c_str())
         ("Solver,S", prog_opt::value(&Solver),
          FormatDefault("Eigensoler to use ("
			+ boost::algorithm::join(LocalEigensolver::EnumerateSolvers(), ", ") + ")", Solver).c_str())
         ("orthogonal", prog_opt::value<std::vector<std::string> >(),
          "force the wavefunction to be orthogonal to this state ***NOT YET IMPLEMENTED***")
         ("dgks", prog_opt::bool_switch(&UseDGKS), "Use DGKS correction for the orthogonality vectors")
         ("oversample", prog_opt::value(&Oversampling.Scale), FormatDefault("For random SVD, oversample by this factor", Oversampling.Scale).c_str())
         ("oversample-min", prog_opt::value(&Oversampling.Add), FormatDefault("For random SVD, minimum amount of oversampling", Oversampling.Add).c_str())
         ("shift-invert-energy", prog_opt::value(&ShiftInvertEnergy), "For the shift-invert and shift-invert-direct solver, the target energy")
         ("subspacesize", prog_opt::value(&SubspaceSize), FormatDefault("Maximum Krylov subspace size for shift-invert solver", SubspaceSize).c_str())
         ("precondition", prog_opt::bool_switch(&UsePreconditioning), "use diagonal preconditioning in the shift-invert solver")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity (can be used more than once)")
         ("flush", prog_opt::bool_switch(&Flush), "Flush terminal output after each line")
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
      BenchFile.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      if (Bench)
      {
         print_preamble(BenchFile, argc, argv);
         BenchFile << "#Time #SweepNum #Site #States #Energy #Trunc #Fidelity #Iter #Tol\n";
      }

      std::cout << "Starting DMRG...\n";
      std::cout << "Input wavefunction: " << FName << std::endl;

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

      dmrg.PreExpansionAlgo = PreExpansionAlgorithm(PreExpandAlgo);
      dmrg.PostExpansionAlgo = PostExpansionAlgorithm(PostExpandAlgo);
      dmrg.UseDGKS = UseDGKS;
      dmrg.Oversampling = Oversampling;
      dmrg.ProjectTwoSiteTangent = ProjectTwoSiteTangent;

      // Defaults for the pre- and post-expansion per sector
      if (!vm.count("pre-sector"))
      {
         if (dmrg.PreExpansionAlgo == PreExpansionAlgorithm::SVD || dmrg.PreExpansionAlgo == PreExpansionAlgorithm::RSVD)
         {
            PreExpand.ExpandPerSector = 0;
         }
         else if (dmrg.PreExpansionAlgo == PreExpansionAlgorithm::RangeFinding || dmrg.PreExpansionAlgo == PreExpansionAlgorithm::Random)
         {
            PreExpand.ExpandPerSector = 1;
         }
      }
      if (!vm.count("post-sector"))
      {
         if (dmrg.PostExpansionAlgo == PostExpansionAlgorithm::SVD || dmrg.PostExpansionAlgo == PostExpansionAlgorithm::RSVD)
         {
            PostExpand.ExpandPerSector = 0;
         }
         else if (dmrg.PostExpansionAlgo == PostExpansionAlgorithm::RangeFinding || dmrg.PostExpansionAlgo == PostExpansionAlgorithm::Random)
         {
            PostExpand.ExpandPerSector = 1;
         }
      }

      dmrg.Solver().SetSolver(Solver);
      dmrg.Solver().MaxTol = MaxTol;
      dmrg.Solver().MinTol = MinTol;
      dmrg.Solver().MinIter = MinIter;
      dmrg.Solver().MaxIter = NumIter;
      dmrg.Solver().FidelityScale = FidelityScale;
      dmrg.Solver().Verbose = Verbose;
      dmrg.Solver().EvolveDelta = EvolveDelta;
      dmrg.Solver().SetShiftInvertEnergy(ShiftInvertEnergy);
      dmrg.Solver().SetSubspaceSize(SubspaceSize);
      dmrg.Solver().SetPreconditioning(UsePreconditioning);
      dmrg.Solver().SetInitialFidelity(1, InitialFidelity);

      EigenSortByWeight = true; // Global variable in density.cpp, to change the eigenvalue sort function

      StatesInfo States;
      States.MinStates = 1;
      States.MaxStates = MaxStates;
      States.TruncationCutoff = TruncCutoff;
      States.EigenvalueCutoff = EigenCutoff;
      std::cout << States << '\n';

      StatesList MyStates(StatesString);
      if (vm.count("steps") && MyStates.size() == 1)
      {
         MyStates.Repeat(NumSweeps);
      }
      std::cout << MyStates;
      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      std::cout << "Number of Lanczos iterations: " << NumIter << std::endl;
      std::cout << "Using solver: " << Solver << std::endl;
      std::cout << "Using pre-expansion algorithm: " << dmrg.PreExpansionAlgo.Name();
      if (dmrg.PreExpansionAlgo != PreExpansionAlgorithm::NoExpansion)
         std::cout << " with " << PreExpand;
      std::cout << "\nUsing post-expansion algorithm: " << dmrg.PostExpansionAlgo.Name();
      if (dmrg.PostExpansionAlgo != PostExpansionAlgorithm::NoExpansion)
         std::cout << " with " << PostExpand;
      std::cout << std::endl;

      int NumStatesKeepNext = MyStates[0].NumStates;
      int ZeroEnvCount = 0;
      double ModFactor = 1.0;
      for (int Sweeps = 0; Sweeps < MyStates.size(); ++Sweeps)
      {
         States.MaxStates = MyStates[Sweeps].NumStates;
         if (MyStates[Sweeps].ZeroEnv)
         {
            ++ZeroEnvCount;
            int ZeroEnvRemain = 0;
            for (int s = Sweeps+1; s < MyStates.size() && MyStates[s].ZeroEnv; ++s)
               ++ZeroEnvRemain;
            ModFactor = double(ZeroEnvRemain) / double(ZeroEnvRemain+ZeroEnvCount);
         }
         NumStatesKeepNext = (!NoGreedy && Sweeps < MyStates.size()-1) ? MyStates[Sweeps+1].NumStates : MyStates[Sweeps].NumStates;
         if (Sweeps % 2 == 0)
            SweepLeft(dmrg, Sweeps+1, States, PreExpand*ModFactor, PostExpand*ModFactor, NumStatesKeepNext);
         else
            SweepRight(dmrg, Sweeps+1, States, PreExpand*ModFactor, PostExpand*ModFactor, NumStatesKeepNext);
      }

      // finished the iterations.
      std::cout << "Orthogonalizing wavefunction...\n";
      Wavefunction.Wavefunction() = dmrg.Wavefunction();

      // any other attributes?
      Wavefunction.Attributes()["LastEnergy"] = formatting::format_complex(dmrg.Solver().LastEnergy());
      Wavefunction.SetDefaultAttributes();

      // History log
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));

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
