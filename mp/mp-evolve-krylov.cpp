// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-evolve-krylov.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
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
#include "common/environment.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int GuessNumKrylov = 8;
      int MinStates = 50;
      int MaxStates = 5000;
      double MixFactor = 0.01;
      bool TwoSite = false;
      std::complex<double> Timestep = 1;
      double EBound = 1e-4;
      double RealTime = 0;
      double Beta = 0;
      int NumTimesteps = 1;
      int SaveTimesteps = 1;
      double ErrorScaleFactor = 0.2;
      bool DoReductionSweep = true;
      std::string OutputPrefix = "";
      bool Normalize = false;
      double KrylovCutoffFactor = 0.2;
      int MaxSweeps = 0;
      int MinKrylovVectors = 3;
      int MaxKrylovVectors = 10;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(), 
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(), "initial wavefunction (required)")
	 ("two-site,2", "modify 2 neighboring sites at once")
	 ("min-states,m", prog_opt::value<int>(&MinStates), "Minimum number of states to keep [default 50]")
	 ("max-states,x", prog_opt::value<int>(&MaxStates), "Maximum number of states to keep [default 5000]")
         ("min-krylov", prog_opt::value(&MinKrylovVectors),
          ("Minimum number of Krylov vectors to use.  The program will abort if fewer than this number of "
           "vectors are relevant - increase the timestep if this happens [default "
           +boost::lexical_cast<std::string>(MinKrylovVectors)+"]").c_str())
         ("max-krylov", prog_opt::value(&MaxKrylovVectors),
          ("Maximum number of Krylov vectors to use.  The program will adaptively reduce the timestep if "
           "this happens [default "
           +boost::lexical_cast<std::string>(MaxKrylovVectors)+"]").c_str())
	 ("mix-factor,f", prog_opt::value<double>(&MixFactor), 
          "Mixing coefficient for the Krylov vector density matrix [default 0.01]")
	 ("timestep,t", prog_opt::value<std::complex<double> >(&Timestep), "Timestep (can be complex) [default 1]")
         ("num-timesteps,n", prog_opt::value<int>(&NumTimesteps), "Total number of timesteps to calculate [default 1]")
         ("save-timesteps,s", prog_opt::value<int>(&SaveTimesteps), 
          "Save the wavefunction after every s timesteps [default 1]")
         ("Time,T", prog_opt::value<double>(&RealTime), 
          "Absolute real time of the input wavefunction (wavefunction attribute \"Time\")")
	 ("Beta,B", prog_opt::value<double>(&Beta),
	  "Absolute imaginary time (inverse temperature) of the input wavefunction "
	  "(wavefunction attribute \"Beta\")")
         ("error-bound,e", prog_opt::value<double>(&EBound), "Truncation error bound per timestep [default 1e-4]")
         ("variance-mod,r", prog_opt::value<double>(&ErrorScaleFactor), 
          "Scale factor for the variance of the Krylov vectors [default 0.2]")
         ("krylov-cutoff", prog_opt::value<double>(&KrylovCutoffFactor),
          "Buffer factor for the cutoff contribution from the Krylov vectors [default 0.2]")
         ("out,o", prog_opt::value<std::string>(&OutputPrefix),
          "Filename prefix for saved wavefunctions [defaults to the initial wavefunction name]")
         ("krylov,k", prog_opt::value<int>(&GuessNumKrylov), "Initial guess for the number of krylov vectors to use [default 8]")
         ("max-sweeps", prog_opt::value(&MaxSweeps), 
          "Maximum number of sweeps to use in the Krylov optimization (dangerous!) [default 0 = no limit]")
         ("full-ortho,F", "full orthogonalization of the Krylov subspace")
         ("normalize", "Normalize the saved wavefunctions.")
         ("reduction-sweep", "Do an additional set of sweeps to try to reduce the dimension of the Krylov vectors")
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

      std::cout << "Starting Krylov...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

      bool FullOrtho = vm.count("full-ortho");
      if (FullOrtho)
         std::cout << "Using full orthogonalization of the Krylov subspace\n";

      DoReductionSweep = vm.count("reduction-sweep");
      if (DoReductionSweep)
         std::cout << "Performing krylov state reduction sweep.\n";

      if (OutputPrefix.empty())
         OutputPrefix = InputWavefunction;

      if (SaveTimesteps != 0)
      {
         std::cout << "Saving the wavefunction every ";
         if (SaveTimesteps == 1)
            std::cout << "timestep.\n";
         else
            std::cout << SaveTimesteps << " timesteps.\n";
         std::cout << "Filename prefix is " << OutputPrefix << '\n';
      }

      Normalize = vm.count("normalize");
      if (Normalize)
         std::cout << "Normalizing the saved wavefunctions.\n";

      // Open the wavefunction
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> P;
      if (SaveTimesteps == 0)
      {
         P = pheap::OpenPersistent(InputWavefunction, CacheSize);
      }
      else
      {
         int NumPageFiles = 1;
         int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
         std::string OutputDir = getenv_or_default("MP_BINPATH", std::string(""));
         if (!OutputDir.empty() && OutputDir[OutputDir.size()-1] != '/')
            OutputDir += '/';

         pheap::Initialize(OutputDir+OutputPrefix+".bin", NumPageFiles, PageSize, CacheSize);
         P = pheap::ImportHeap(InputWavefunction);
      }
      MPWavefunction Psi = *P;
      P = pvalue_ptr<MPWavefunction>();
      // make sure the wavefunction is normalized ** we shouldn't need to do this? **
      // Psi.normalize();

#if 0
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
#endif

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
      KrylovLoop dmrg((SplitOperator(Hamiltonian)), 
                      SplitOperator(prod(Hamiltonian, Hamiltonian, Hamiltonian.TransformsAs())), 
                      CenterWavefunction(Psi));

      TwoSite = vm.count("two-site");
      if (TwoSite)
      {
	 std::cout << "Optimizing two sites at a time" << std::endl;
      }

      std::cout << "Density matrix mixing coefficient: " << MixFactor << std::endl;
      //      std::cout << "Number of half-sweeps: " << NumSweeps << std::endl;
      std::cout << "Timestep: " << Timestep << std::endl;
      std::cout << "Error bound per timestep: " << EBound << std::endl;

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
      dmrg.ErrorScaleFactor_ = ErrorScaleFactor;
      dmrg.LastNumKrylov_ = GuessNumKrylov;
      dmrg.DoReductionSweep_ = DoReductionSweep;
      dmrg.KrylovCutoffFactor_ = KrylovCutoffFactor;
      dmrg.MaxSweeps_ = (MaxSweeps > 0 ? MaxSweeps : 1000000);

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = MaxStates;
      SInfo.TruncationCutoff = EBound / Psi.size(); // truncation per lattice site
      //SInfo.EigenvalueCutoff = EigenCutoff;
            
      for (int CurrentTimestep = 0; CurrentTimestep < NumTimesteps; ++CurrentTimestep)
      {
         dmrg.ConstructKrylovBasis(Timestep);
         std::cout << "Constructing evolved wavefunction..." << std::endl;
         dmrg.Evolve(Timestep,SInfo);
         std::cout << "done.  Current time is " << dmrg.RealTime_
                   << ", current beta is " << dmrg.Beta_ << std::endl;
         TRACE(dmrg.Norm());

         if (dmrg.LastNumKrylov_ < MinKrylovVectors)
         {
            std::cerr << "fatal: number of Krylov vectors (" << dmrg.LastNumKrylov_ 
                      << ") has dropped below the minimum possible (" << MinKrylovVectors << ").  Aborting.\n";
            std::exit(2);
         }

         if (SaveTimesteps > 0 && (CurrentTimestep+1) % SaveTimesteps == 0)
         {
            std::ostringstream filetime;
            filetime.precision(5);
            if (dmrg.RealTime_ != 0.)
               filetime << ".t" << dmrg.RealTime_;
            if (dmrg.Beta_ != 0.)
               filetime << ".b" << dmrg.Beta_;
            std::string FName = OutputPrefix + filetime.str();
            std::cout << "Saving wavefunction as " << FName << '\n';
            pvalue_ptr<MPWavefunction> OutPsi(new MPWavefunction(dmrg.Wavefunction().AsLinearWavefunction()));
            if (Normalize)
               OutPsi.mutate()->normalize();
            pheap::ExportHeap(FName, OutPsi);
         }

         if (dmrg.LastNumKrylov_ > MaxKrylovVectors)
         {
            Timestep *= 0.5;
            CurrentTimestep *= 2;
            NumTimesteps *= 2;
            SaveTimesteps *= 2;
            ++CurrentTimestep; // Since we just halved Timestep, we have effectively done two timesteps this iteration
            std::cout << "Number of Krylov vectors (" << dmrg.LastNumKrylov_ 
                      << ") has exceeded the maximum (" << MaxKrylovVectors << ".  Decreasing timestep to "
                      << Timestep << std::endl;
         }
      }

      if (SaveTimesteps == 0)
      {
         P = pvalue_ptr<MPWavefunction>(new MPWavefunction(dmrg.Wavefunction().AsLinearWavefunction()));
         if (Normalize)
            P.mutate()->normalize();
         pheap::ShutdownPersistent(P);
      }
      else
         pheap::Shutdown();

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
