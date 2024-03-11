// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-evolve-magnus.cpp
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

using std::sin;

// global parameters for the Hamiltonian; these are set from the command line
double Period;     // period of the hopping oscillation
double Delta;      // amplitude of the hopping oscillation
double U;          // coulomb term
std::string Hopping = "Hopping";  // default operator name for the hopping Hamiltonain
std::string Coulomb = "Coulomb";  // default operator name for the interaction

MPOperator Hamiltonian(OperatorList const& Lattice, double Time)
{
   TRACE(Time)(Period)(Delta)(U);
   return (1.0 + Delta * sin(2.0 * math_const::pi * Time / Period)) * Lattice[Hopping]
      + U * Lattice[Coulomb];
}

MPOperator Omega(OperatorList const& Lattice, double Time, double Timestep, int Order)
{
   if (Order == 2)
   {
      // Result = H(t + 0.5 \delta t)
      return Timestep * Hamiltonian(Lattice, Time + 0.5*Timestep);
   }
   else if (Order == 4)
   {
      // H1 = H(t + delta t * (0.5 - sqrt(3)/6))
      MPOperator H1 = Hamiltonian(Lattice, Time + Timestep * (0.5 - std::sqrt(3.0)/6.0));
      // H2 = H(t + delta t * (0.5 + sqrt(3)/6))
      MPOperator H2 = Hamiltonian(Lattice, Time + Timestep * (0.5 + std::sqrt(3.0)/6.0));

      // Commutator = [H2, H1]
      MPOperator Commutator = prod(H2, H1, H1.TransformsAs())
         - prod(H1, H2, H1.TransformsAs());

      // Result = (Timestep/2) * (H1 + H2) + (Timestep*sqrt(3)/12) * [H2, H1]
      return 0.5*Timestep*(H1 + H2) +
         std::complex<double>(0.0, Timestep*Timestep*std::sqrt(3.0)/12.0) * Commutator;
   }
   else
   {
      PANIC("fatal: invalid Magnus expansion order, must be 2 or 4")(Order);
      return MPOperator();
   }
}

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int MinStates = 50;
      int MaxStates = 5000;
      double MixFactor = 0.01;
      bool TwoSite = false;
      double Timestep = 0;
      double Time = 0;
      double RealTime = 0;
      double EBound = 1e-4;
      int Order = 4;
      std::string LatticeFile;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("lattice,l", prog_opt::value<std::string>(&LatticeFile),
          "Lattice file to use to extract the Hamiltonian")
         ("wavefunction,w", prog_opt::value<std::string>(),
          "initial wavefunction (required)")
         ("U", prog_opt::value<double>(&U),
          "Coulomb interaction (required)")
         ("delta,D", prog_opt::value<double>(&Delta),
          "Amplitude of the hopping oscillation (required)")
         ("period,p", prog_opt::value<double>(&Period),
          "period of oscillation of the hopping term")
         ("step,s", prog_opt::value<double>(&Timestep), "Timestep (required)")
         ("time,t", prog_opt::value<double>(&Time),
          "Time to evolve to (absolute value)")
         ("numsteps,n", prog_opt::value<int>(),
          "Number of timesteps to evolve (alternative to --time)")
         ("order", prog_opt::value<int>(&Order),
          "Order of the Magnus expansion to use, valid values are 2 or 4 [default 4]")
         ("two-site,2", "modify 2 neighboring sites at once")
         ("min-states,m", prog_opt::value<int>(&MinStates),
          "Minimum number of states to keep [default 50]")
         ("max-states,x", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 5000]")
         ("mix-factor,f", prog_opt::value<double>(&MixFactor),
          "Mixing coefficient for the density matrix [default 0.01]")
         ("Time,T", prog_opt::value<double>(&RealTime),
          "Absolute real time of the input wavefunction "
          "(wavefunction attribute \"Time\", default 0)")
         ("error-bound,e", prog_opt::value<double>(&EBound),
          "Error bound per unit time [default 1e-4]")
         ("full-ortho,F", "full orthogonalization of the Krylov subspace")
         ("hopping", prog_opt::value<std::string>(&Hopping),
          "Operator name for the hopping operator")
         ("coulomb", prog_opt::value<std::string>(&Coulomb),
          "Operator name for the on-site Coulomb operator")
          ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || !vm.count("wavefunction")
          || !vm.count("U") || !vm.count("delta") || !vm.count("period")
          || !vm.count("step") || !vm.count("lattice"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-evolve-magnus [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (vm.count("numsteps"))
      {
         if (vm.count("time"))
         {
            std::cerr << "fatal: only one option out of \"--numsteps\" "
               "and \"--time\" may be specified!\n";
            return 1;
         }
         Time = Timestep * vm["numsteps"].as<int>();
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

      // Get the initial time from the attributes, only if it is not specified on the command line
      if (vm.count("Time") == 0 && Psi.Attributes().count("Time"))
      {
         RealTime = Psi.Attributes()["Time"].as<double>();
      }

      // Set up the lattice and Omega operator
      pvalue_ptr<OperatorList> Lattice = pheap::ImportHeap(LatticeFile);
      MPOperator Ham = Omega(*Lattice, RealTime, Timestep, Order);

      TRACE(expectation(Psi, Ham, Psi));
      TRACE(expectation(Psi, prod(Ham, Ham, Ham.TransformsAs()), Psi));
      // Now we can construct the actual KrylovSolver object
      KrylovLoop dmrg(Ham, prod(Ham, Ham, Ham.TransformsAs()), Psi);

      TwoSite = vm.count("two-site");
      if (TwoSite)
      {
         std::cout << "Optimizing two sites at a time" << std::endl;
      }

      std::cout << "Density matrix mixing coefficient: " << MixFactor << std::endl;
      //      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      std::cout << "Timestep: " << Timestep << std::endl;
      std::cout << "Error bound per unit time: " << EBound << std::endl;

      // The timestep is measured in units of -i/hbar
      // get the normalized timestep
      double TimeDistance = LinearAlgebra::norm_2(Timestep);
      std::complex<double> TimeDir = Timestep * std::complex<double>(0.0,-1.0) / TimeDistance;
      dmrg.TimeDirection_ = TimeDir;
      dmrg.ErrorBound_ = EBound;
      dmrg.MinStates_ = MinStates;
      dmrg.MaxStates_ = MaxStates;
      dmrg.FullOrtho_ = FullOrtho;
      dmrg.MixFactor_ = MixFactor;
      dmrg.TwoSite_ = TwoSite;
      dmrg.RealTime_ = RealTime;
      dmrg.Beta_ = 0;

      // Get a lower bound on the number of steps to do
      int NumSteps = int(floor(Time / TimeDistance
                           - 100*std::numeric_limits<double>::epsilon()));
      TRACE(NumSteps);
      if (NumSteps < 0)
         NumSteps = 0;

      for (int i = 0; i < NumSteps; ++i)
      {
         dmrg.ConstructKrylovBasis(Timestep);
         std::cout << "Constructing evolved wavefunction..." << std::endl;
         dmrg.Evolve(Timestep, EBound * Timestep,
                     Omega(*Lattice, dmrg.RealTime_ + Timestep, Timestep, Order));
         std::cout << "done.  Current time is " << dmrg.RealTime_ << std::endl;
         TRACE(dmrg.Norm());
      }

      // and the final timestep
      double FinalTimestep = Time + RealTime - dmrg.RealTime_;
      TRACE(FinalTimestep);
      dmrg.ConstructKrylovBasis(FinalTimestep);
      std::cout << "Constructing evolved wavefunction..." << std::endl;
      dmrg.Evolve(FinalTimestep, EBound * Timestep);
      dmrg.RealTime_ = RealTime+Time;  // to avoid messy rounding issues, just
      // set the time to what we expect it to be
      std::cout << "done.  Current time is " << dmrg.RealTime_ << std::endl;
      TRACE(dmrg.Norm());

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
