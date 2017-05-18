// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-cv2.cpp
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/functional-solver.h"
#include "mp/copyright.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "interface/operator-parser.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include <boost/program_options.hpp>
#include "common/prog_opt_accum.h"
#include "interface/inittemp.h"
#include <iostream>
#include "common/environment.h"

namespace prog_opt = boost::program_options;

// helper function to load an attribute from the command line or wavefunction
template <typename T>
void LoadAttribute(prog_opt::variables_map& Options, MPWavefunction& Psi,
                   std::string const& Name, T& Value, bool Show = false)
{
   if (Options.count(Name) == 1)
      {
         Value = Options[Name]. template as<T>();
         Psi.Attributes()[Name] = Value;
      }
   else
      {
         if (Psi.Attributes().count(Name) == 0)
         {
            std::ostringstream ostr;
            ostr << "fatal: missing parameter: " << Name << ".";
            throw std::runtime_error(ostr.str());
         }
         Value = Psi.Attributes()[Name].template as<T>();
      }
   if (Show)
      std::cout << Name << ": " << Value << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      int NumIter = 20;
      int MinIterations = 0;
      double GroundstateEnergy = 0;
      double Frequency = 0;
      double Broadening = 0;
      double MixFactor = 0.0;
      int NumSweeps = 2;
      int MinStates = 1;
      int MaxStates = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = -1;
      int Verbose = 0;
      std::string LanczosStr, WavefunctionStr, HamStr;
      bool TwoSite = false;
      double Precision = 1E-15;
      bool TwoStepKrylov = false;
      bool MixNormalize = false;
      double LanczosMixFactor = 0.0;
      bool LinearSolver = false;
      bool SquareMeHarder = false;
      bool UseResid = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("lanczos,l", prog_opt::value(&LanczosStr),
          "input Lanczos vector (required)")
         ("wavefunction,w", prog_opt::value(&WavefunctionStr),
          "correction vector (required)")
         ("Frequency,F", prog_opt::value(&Broadening),
          "Frequency (wavefunction attribute \"Frequency\")")
         ("Broadening,B", prog_opt::value(&Broadening),
          "Broadening (wavefunction attribute \"Broadening\")")
         ("GroundstateEnergy,G", prog_opt::value(&GroundstateEnergy),
          "groundstate energy of the Hamiltonian (wavefunction attribute"
          " \"GroundstateEnergy\"")
         ("two-site,2", prog_opt::bool_switch(&TwoSite),
          "modify 2 neighboring sites at once")
         ("iter,i", prog_opt::value(&NumIter),
          ("Maximum size of the Krylov subspace [default "
           + boost::lexical_cast<std::string>(NumIter)+"]").c_str())
         ("min-iter", prog_opt::value(&MinIterations),
          ("Minimum size of the Krylov subspace [default "
           + boost::lexical_cast<std::string>(MinIterations)+"]").c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
          ("Minimum number of states to keep [default "
           +boost::lexical_cast<std::string>(MinStates)+"]").c_str())
         ("max-states,m", prog_opt::value<int>(&MaxStates),
          ("Maximum number of states to keep [default "
           +boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("trunc,r", prog_opt::value(&TruncCutoff),
          ("Cutoff truncation error per site [default "
           +boost::lexical_cast<std::string>(TruncCutoff)+"]").c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          ("Cutoff threshold for density matrix eigenvalues [default "
           +boost::lexical_cast<std::string>(EigenCutoff)+"]").c_str())
         ("mix-factor,f", prog_opt::value(&MixFactor),
          ("Mixing coefficient for the density matrix [default "
           + boost::lexical_cast<std::string>(MixFactor)+"]").c_str())
         ("precision,p", prog_opt::value(&Precision),
          ("Precision for the functional solver - stop if the change in one "
           "iteration is less than this number [default "
           + boost::lexical_cast<std::string>(Precision)+"]").c_str())
         ("twostepkrylov", prog_opt::bool_switch(&TwoStepKrylov),
          "Add krylov vectors H|k_n> and H^2|k_{n-1}> alternately")
         ("lanczos-mix", prog_opt::value(&LanczosMixFactor),
          ("Mix the Lanczos vector into the basis with this weight [default "
            + boost::lexical_cast<std::string>(LanczosMixFactor)+"]").c_str())
         ("mix-normalize", prog_opt::bool_switch(&MixNormalize),
          "Normalize the density matrices for the cv and the lv to 1 before mixing")
         ("linear-solver", prog_opt::bool_switch(&LinearSolver),
          "Use a linear solver instead of the functional minimization")
         ("use-resid", prog_opt::bool_switch(&UseResid),
          "Construct the Krylov basis of the functional minimiziation using the best approximate residual vector")
         ("square", prog_opt::bool_switch(&SquareMeHarder),
          "Don't use the exact matrix elements of H^2 in the functional minimization, instead calculate H twice")
         ("sweeps,s", prog_opt::value(&NumSweeps),
          ("Number of half-sweeps to perform [default "
           + boost::lexical_cast<std::string>(NumSweeps)+"]").c_str())
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "increase verbosity (can be used more than once)")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("lanczos") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-cv2 [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(WavefunctionStr, 655360);
      MPWavefunction Psi = *PsiPtr;
      PsiPtr = pvalue_ptr<MPWavefunction>();

      // Set up the Hamiltonian
      LoadAttribute(vm, Psi, "Hamiltonian", HamStr);
      if (Verbose >= 1)
         std::cerr << "Hamiltonian: " << HamStr << std::endl;
      OperatorList OpList;
      MPOperator H;
      std::tie(OpList, H) = ParseLatticeAndOperator(HamStr);

      LoadAttribute(vm, Psi, "GroundstateEnergy", GroundstateEnergy, Verbose >= 1);
      LoadAttribute(vm, Psi, "Frequency", Frequency, Verbose >= 1);
      LoadAttribute(vm, Psi, "Broadening", Broadening, Verbose >= 1);

      if (TwoSite && Verbose >= 1)
         std::cout << "Modifying two sites per step.\n";

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;
      if (Verbose >= 1)
         std::cout << SInfo << '\n';

      pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(LanczosStr);

      MPOperator Part = (Frequency+GroundstateEnergy) * OpList["I"] - H;
      MPOperator A = Part*Part + Broadening*Broadening*OpList["I"];
      MPOperator B = Part + std::complex<double>(0.0, Broadening)*OpList["I"];
      MPOperator Z = -1.0*(Part + std::complex<double>(0.0, -Broadening)*OpList["I"]);
      MPOperator Zb = -1.0*(Part + std::complex<double>(0.0, Broadening)*OpList["I"]);

      double Shift = norm_frob_sq(*Rhs);

      FunctionalSolver solver(CenterWavefunction(Psi), A, B, Z, Zb,
                              CenterWavefunction(*Rhs), Frequency, Broadening);
      Psi = MPWavefunction();
      solver.Precision = Precision;
      solver.TwoStepKrylov = TwoStepKrylov;
      solver.MixNormalize = MixNormalize;
      solver.LanczosMixFactor = LanczosMixFactor;
      solver.MinIterations = MinIterations;
      solver.LinearSolver = LinearSolver;
      solver.SquareMeHarder = SquareMeHarder;
      solver.UseResid = UseResid;

      bool First = false;
      for (int Sweeps = 0; Sweeps < NumSweeps/2; ++Sweeps)
      {

         {
            solver.ExpandLeft();
            double E = First ? 0.0 : solver.Solve(NumIter) + Shift;
            TruncationInfo States = solver.TruncateLeft(SInfo, MixFactor);

            std::cout << "Partition=(" << solver.LeftSize() << ',' << solver.RightSize()
                      << ") r2=" << E
                      << " States=" << States.KeptStates()
                      << " Trunc=" << States.TruncationError()
                      << " NIter=" << solver.IterationNumMultiplies
                      << '\n';
         }

         // sweep right
         while (solver.RightSize() > 1)
         {
            solver.ShiftRightAndExpand();
            if (TwoSite)
               solver.ExpandRight();
            double E = First ? 0.0 : solver.Solve(NumIter) + Shift;
            TruncationInfo States = solver.TruncateLeft(SInfo, MixFactor);

            std::cout << "Partition=(" << solver.LeftSize() << ',' << solver.RightSize()
                      << ") r2=" << E
                      << " States=" << States.KeptStates()
                      << " Trunc=" << States.TruncationError()
                      << " NIter=" << solver.IterationNumMultiplies
                      << '\n';
         }
         First = false;


   //   std::cout << solver.ResidualNorm() << ' ' << solver.ExactResidualNorm() << '\n';
   //std::cout //<< (difference(solver.Wavefunction(), xOld) / norm_2(xOld)) << ' '
      //             << (difference(solver.WavefunctionAx(), AxOld) / norm_2(AxOld)) << ' '
      //           << (overlap(solver.Wavefunction(), xOld) / norm_2_sq(xOld)) << '\n';
   //xOld = solver.Wavefunction();
   //   AxOld = solver.WavefunctionAx();

         {
            solver.ExpandRight();
            double E = solver.Solve(NumIter) + Shift;
            TruncationInfo States = solver.TruncateRight(SInfo, MixFactor);

            std::cout << "Partition=(" << solver.LeftSize() << ',' << solver.RightSize()
                      << ") r2=" << E
                      << " States=" << States.KeptStates()
                      << " Trunc=" << States.TruncationError()
                      << " NIter=" << solver.IterationNumMultiplies
                      << '\n';
         }

         // sweep left
         while (solver.LeftSize() > 1)
         {
            solver.ShiftLeftAndExpand();
            if (TwoSite)
               solver.ExpandLeft();
            double E = solver.Solve(NumIter) + Shift;
            TruncationInfo States = solver.TruncateRight(SInfo, MixFactor);

            std::cout << "Partition=(" << solver.LeftSize() << ',' << solver.RightSize()
                      << ") r2=" << E
                      << " States=" << States.KeptStates()
                      << " Trunc=" << States.TruncationError()
                      << " NIter=" << solver.IterationNumMultiplies
                      << '\n';
         }

   //   std::cout << solver.ResidualNorm() << ' ' << solver.ExactResidualNorm() << '\n';
   // std::cout //<< (difference(solver.Wavefunction(), xOld) / norm_2(xOld)) << ' '
      //             << (difference(solver.WavefunctionAx(), AxOld) / norm_2(AxOld)) << ' '
      //          << (overlap(solver.Wavefunction(), xOld) / norm_2_sq(xOld)) << '\n';
   // xOld = solver.Wavefunction();
   //   AxOld = solver.WavefunctionAx();

      }

      PsiPtr = new MPWavefunction(solver.Wavefunction().AsLinearWavefunction());
      pheap::ShutdownPersistent(PsiPtr);

   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
