// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-evolve-krylov-x.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/krylovloop.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "linearalgebra/matrix_utility.h"
#include <iostream>
#include <cmath>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int NumKrylov = 4;
      int MinStates = 50;
      int MaxStates = 5000;
      double MixFactor = 0.01;
      bool TwoSite = false;
      double HNorm = -1;
      double InitialTimestep = 0.25;
      std::complex<double> Timestep = 1;
      double EBound = 1e-4;
      double RealTime = 0;
      double Beta = 0;
      double MaxTimestep = 1e10;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(),
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "initial wavefunction (required)")
         ("two-site,2", "modify 2 neighboring sites at once")
         ("krylov,k", prog_opt::value<int>(&NumKrylov), "Number of Krylov vectors [default 4]")
         ("min-states,m", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 50]")
         ("max-states,x", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 5000]")
         ("mix-factor,f", prog_opt::value<double>(&MixFactor),
          "Mixing coefficient for the density matrix [default 0.01]")
         ("hnorm,h", prog_opt::value<double>(&HNorm),
          "2-norm of the Hamiltonian operator (used to guess the timestep for the first iteration)")
         ("initialstep,i", prog_opt::value<double>(&InitialTimestep),
          "Initial timestep to use (alternative to specifying --hnorm)")
         ("time,t", prog_opt::value<std::complex<double> >(&Timestep), "Time to evolve (can be complex) [default 1]")
         ("Time,T", prog_opt::value<double>(&RealTime),
          "Absolute real time of the input wavefunction (wavefunction attribute \"Time\")")
         ("Beta,B", prog_opt::value<double>(&Beta),
          "Absolute imaginary time (inverse temperature) of the input wavefunction "
          "(wavefunction attribute \"Beta\")")
         ("error-bound,e", prog_opt::value<double>(&EBound), "Error bound per unit time [default 1e-4]")
         ("max-timestep", prog_opt::value<double>(&MaxTimestep),
          "Maximum timestep [default 1e10]")
         ("full-ortho,F", "full orthogonalization of the Krylov subspace")
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
         std::cerr << "usage: mp-evolve-krylov [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      bool FullOrtho = vm.count("full-ortho");


      std::cout << "Starting Krylov...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      // Open the wavefunction
      pvalue_ptr<MPWavefunction> P = pheap::OpenPersistent(InputWavefunction, 655360);
      MPWavefunction Psi = *P;
      P = pvalue_ptr<MPWavefunction>();
      // make sure the wavefunction is normalized
      Psi.normalize();

      // Make sure the center matrix is at one edge
      if (Psi.LeftSize() != 1 && Psi.RightSize() != 1)
      {
         TRACE(Psi.LeftSize())(Psi.RightSize());
         std::cout << "The center matrix is not located at an edge.  Rotating..." << std::flush;
         if (Psi.LeftSize() > Psi.RightSize())
         {
            while (Psi.RightSize() > 1)
               Psi.RotateRight();
         }
         else
         {
            while (Psi.LeftSize() > 1)
               Psi.RotateLeft();
         }
         std::cout << "done" << std::endl;
      }

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
         Psi.Attributes()["Hamiltonian"] = HamString;
      }
      else
      {
         if (Psi.Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the initial wavefunction.\n";
            return 1;
         }
         HamString = Psi.Attributes()["Hamiltonian"].as<std::string>();
      }
      std::cout << "Hamiltonian: " << HamString << std::endl;
      MPOperator Hamiltonian = ParseOperator(HamString);

      // Get the initial time from the attributes, only if it is not specified on the command line
      if (vm.count("Time") == 0 && Psi.Attributes().count("Time"))
      {
         RealTime = Psi.Attributes()["Time"].as<double>();
      }
      // Get the initial beta from the attributes, only if it is not specified on the command line
      if (vm.count("Beta") == 0 && Psi.Attributes().count("Beta"))
      {
         Beta = Psi.Attributes()["Beta"].as<double>();
      }

      // Now we can construct the actual KrylovSolver object
      KrylovLoop dmrg(Hamiltonian, prod(Hamiltonian, Hamiltonian, Hamiltonian.TransformsAs()), Psi);

      TwoSite = vm.count("two-site");
      if (TwoSite)
      {
         std::cout << "Optimizing two sites at a time" << std::endl;
      }

      std::cout << "Density matrix mixing coefficient: " << MixFactor << std::endl;
      std::cout << "Number of Krylov vectors: " << NumKrylov << std::endl;
      //      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      std::cout << "Timestep: " << Timestep << std::endl;
      std::cout << "Error bound per unit time: " << EBound << std::endl;

      // The timestep is measured in units of -i/hbar
      Timestep *= std::complex<double>(0.0, -1.0);

      // get the normalized timestep
      double TimeDistance = LinearAlgebra::norm_2(Timestep);
      std::complex<double> TimeDir = Timestep / TimeDistance;
      dmrg.TimeDirection_ = TimeDir;
      dmrg.ErrorBound_ = EBound;
      dmrg.MinStates_ = MinStates;
      dmrg.MaxStates_ = MaxStates;
      dmrg.FullOrtho_ = FullOrtho;
      dmrg.MixFactor_ = MixFactor;
      dmrg.TwoSite_ = TwoSite;
      dmrg.RealTime_ = RealTime;
      dmrg.Beta_ = Beta;

      double CurrentTime = 0;

      // Determine the first timestep, either specified as a parameter or
      // from the Hochbruck bound given the 2-norm of the Hamiltonian
      double LastTimestep = InitialTimestep;
      if (HNorm > 0)
         LastTimestep = dmrg.TimestepFromHochbruckBound(NumKrylov, HNorm, std::sqrt(EBound));
      std::cout << "Guessed initial timestep is " << LastTimestep << '\n';
      TRACE(LastTimestep);
      while (CurrentTime < TimeDistance)
      {
         if (LastTimestep > TimeDistance-CurrentTime)
            LastTimestep = TimeDistance-CurrentTime;
         dmrg.ConstructKrylovBasis(NumKrylov, LastTimestep);
         double NextTimestep = dmrg.GetTimestep(LastTimestep);
         if (NextTimestep > MaxTimestep)
            NextTimestep = MaxTimestep;

         // Dont evolve beyond where we want to go!
         if (NextTimestep >= TimeDistance - CurrentTime)
         {
            NextTimestep = TimeDistance - CurrentTime;
            CurrentTime = TimeDistance;  // this should ensure the loop terminates immediately
         }
         else
         {
            CurrentTime += NextTimestep;
         }
         std::cout << "NextTimestep = " << NextTimestep << '\n';
         std::cout << "Constructing evolved wavefunction..." << std::endl;
         dmrg.Evolve(NextTimestep, EBound * EBound * NextTimestep * NextTimestep);
         std::cout << "done.  Current time is " << dmrg.RealTime_
                   << ", current beta is " << dmrg.Beta_ << std::endl;
         LastTimestep = NextTimestep;
         TRACE(dmrg.Norm());
      }

      P = pvalue_ptr<MPWavefunction>(new MPWavefunction(dmrg.Wavefunction()));
      pheap::ShutdownPersistent(P);

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
