// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-evolve-krylov-simple.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/splitoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/simplekrylov.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "linearalgebra/matrix_utility.h"
#include <iostream>
#include <cmath>
#include "common/environment.h"
#include <boost/none.hpp>
#include <boost/optional.hpp>

namespace prog_opt = boost::program_options;

typedef std::complex<double> complex;

class KrylovLoop
{
   public:
      KrylovLoop(CenterWavefunction const& Psi, 
                 SplitOperator const& Hamiltonian,
                 complex Timestep, bool UseK1 = false);

      void SetMinStates(int MinStates_) { States.MinStates = MinStates_; }
      void SetMaxStates(int MaxStates_) { States.MaxStates = MaxStates_; }
      void SetTwoSite(bool TwoSite_) { TwoSite = TwoSite_; }
      void SetMinTrunc(double MinTrunc_) { States.TruncationCutoff = MinTrunc_; }
      void SetMaxIter(int MaxIter_) { MaxIter = MaxIter_; }
      void SetConvergenceFactor(double x_) { ConvergenceFactor = x_; }
      void SetMinSweeps(int s) { MinSweeps = s; }
      void SetMaxSweeps(int s) { MaxSweeps = s; }
      void SetMixFactor(double x) { MixFactor = x; }
      void SetMaxEigen(double e) { States.EigenvalueCutoff = e; } 
      void SetSolverTol(double x) { ErrBound = x; }
      void SetUseFinalSweep(bool x) { UseFinal = x; }
      void SetupLogFiles(std::string const& Prefix) { dmrg.SetupLogFiles(Prefix); }

      bool IsConverged() const;

      // Advance the timestep 
      void NextTimestep();

      CenterWavefunction const& Wavefunction() const { return dmrg.Wavefunction(); }

      void DoSweepUntilConverged();

      double Difference() const;

      void DoSweep();

   private:
      void SweepRight();
      void SweepLeft();

      // the number of sweeps we have done this timestep
      int NumSweepsThisTimestep;

      SimpleKrylov dmrg;

      CenterWavefunction PsiPrev;  // the previous wavefunction, for convergence check

      mutable boost::optional<double> Diff; // cache of |Psi-PsiPrev|^2

      int SystemSize;    // number of sites

      int MinSweeps;
      int MaxSweeps;
      double ConvergenceFactor;  // stop when |Psi-PsiPrev|^2 < ConvergenceFactor * Trunc
      int MaxIter;
      bool TwoSite;
      double MixFactor;
      StatesInfo States;
      bool UseFinal;
      double ErrBound;
};

KrylovLoop::KrylovLoop(CenterWavefunction const& Psi, 
                       SplitOperator const& Hamiltonian,
                       complex Timestep, bool UseK1)
   : dmrg(Psi, Hamiltonian, Timestep, Psi, UseK1),
     PsiPrev(Psi),
     SystemSize(Psi.size()),
     MinSweeps(0),
     MaxSweeps(10000),
     ConvergenceFactor(1.0),
     MaxIter(10),
     TwoSite(false),
     MixFactor(0),
     UseFinal(false),
     ErrBound(0.0)
{
}

void KrylovLoop::SweepRight()
{
   dmrg.StartSweep();
   dmrg.StartIteration();
   dmrg.ExpandLeft();
   if (TwoSite) dmrg.ExpandRight();
   double Iter = dmrg.Solve(MaxIter, ErrBound);
   TruncationInfo TInfo = dmrg.TruncateLeft(States, MixFactor);
   dmrg.EndIteration();
   std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ')'
	     << " States:" << TInfo.KeptStates() 
             << " Variance:" << dmrg.IterationOverallVariance
      //<< " Truncation-error:" << TInfo.TruncationError() 
             << " Truncation:" << dmrg.IterationRealTruncation
             << " Iterations:" << Iter << '\n';
   // sweep right
   while (dmrg.RightSize() > 1)
   {
      dmrg.StartIteration();
      dmrg.ShiftRightAndExpand();
      if (TwoSite) dmrg.ExpandRight();
      Iter = dmrg.Solve(MaxIter, ErrBound);
      TInfo = dmrg.TruncateLeft(States, MixFactor);
      dmrg.EndIteration();
      std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ')'
		<< " States:" << TInfo.KeptStates() 
                << " Variance:" << dmrg.IterationOverallVariance
         //<< " Truncation-error:" << TInfo.TruncationError() 
                << " Truncation:" << dmrg.IterationRealTruncation
                << " Iterations:" << Iter << '\n';
   }
   dmrg.EndSweep();
   ++NumSweepsThisTimestep;
}

void KrylovLoop::SweepLeft()
{
   dmrg.StartSweep();
   dmrg.StartIteration();
   dmrg.ExpandRight();
   if (TwoSite) dmrg.ExpandLeft();
   double Iter = dmrg.Solve(MaxIter, ErrBound);
   TruncationInfo TInfo = dmrg.TruncateRight(States, MixFactor);
   dmrg.EndIteration();
   std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ')'
	     << " States:" << TInfo.KeptStates() 
             << " Variance:" << dmrg.IterationOverallVariance
      //<< " Truncation-error:" << TInfo.TruncationError() 
             << " Truncation:" << dmrg.IterationRealTruncation
             << " Iterations:" << Iter << '\n';
   // sweep left
   while (dmrg.LeftSize() > 1)
   {
      dmrg.StartIteration();
      dmrg.ShiftLeftAndExpand();
      if (TwoSite) dmrg.ExpandRight();
      Iter = dmrg.Solve(MaxIter, ErrBound);
      TInfo = dmrg.TruncateRight(States, MixFactor);
      dmrg.EndIteration();
      std::cout << "Partition:(" << dmrg.LeftSize() << ',' << dmrg.RightSize() << ')'
		<< " States:" << TInfo.KeptStates() 
                << " Variance:" << dmrg.IterationOverallVariance
         //<< " Truncation-error:" << TInfo.TruncationError() 
                << " Truncation:" << dmrg.IterationRealTruncation
                << " Iterations:" << Iter << '\n';
   }
   dmrg.EndSweep();
   ++NumSweepsThisTimestep;
}

void KrylovLoop::DoSweep()
{
   Diff = boost::none;
   PsiPrev = dmrg.Wavefunction();
   PsiPrev.normalize();
   if (PsiPrev.LeftSize() == 1)
      this->SweepRight();
   else
      this->SweepLeft();
}

double KrylovLoop::Difference() const
{
   if (!Diff)
   {
      CenterWavefunction w = dmrg.Wavefunction();
      w.normalize();
      Diff = scalar_difference_sq(w.AsLinearWavefunction(), 
                                  PsiPrev.AsLinearWavefunction());
   }
   return Diff.get();
}

bool KrylovLoop::IsConverged() const
{
   if (NumSweepsThisTimestep == 0)
      return false;

   // we could also used the requested truncation instead
   return this->Difference() < dmrg.SweepTruncation * ConvergenceFactor;
}

void KrylovLoop::DoSweepUntilConverged()
{
   while ((!this->IsConverged() || NumSweepsThisTimestep < MinSweeps)
          && NumSweepsThisTimestep < MaxSweeps)
   {
      this->DoSweep();
      std::cout << "Sweep finished, cumulative truncation = " << dmrg.SweepTruncation
                << ", wavefunction difference = " << this->Difference() << std::endl;
      std::cout << "Overall variance = " << dmrg.SweepOverallVariance << "\n";
      std::cout << "Variance bound = " << (dmrg.SweepStDevBound * dmrg.SweepStDevBound) << "\n";
      std::cout << "Sweep real truncation = " << dmrg.SweepRealTruncation << "\n";
   }
   if (this->IsConverged())
      std::cout << "Convergence criteria satisfied.";
   else
      std::cout << "Maximum number of sweeps reached (" << NumSweepsThisTimestep << ")";
   std::cout << std::endl;

   // Final sweep option
   if (UseFinal)
   {
      std::cout << "Doing final sweep with zero mix factor\n";
      StatesInfo OldStates = States;
      if (!TwoSite)
      {
         // zero truncation error doesn't make much sense for TwoSite
         States.TruncationCutoff = 0;
         States.EigenvalueCutoff = 0;
      }
      double OldMixFactor = MixFactor;
      MixFactor = 0;
      this->DoSweep();
      MixFactor = OldMixFactor;
      States = OldStates;
   }
}

void KrylovLoop::NextTimestep()
{
   PsiPrev = dmrg.Wavefunction();
   NumSweepsThisTimestep = 0;
   dmrg.AdvanceTime();
   //   dmrg = SimpleKrylov(PsiPrev, dmrg.Ham, dmrg.Timestep, dmrg.Psi1_H_Psi1);
}

int main(int argc, char** argv)
{
   try
   {
      int MinStates = 10;
      int MaxStates = DefaultMaxStates;
      bool TwoSite = false;
      std::complex<double> Timestep = 1;
      double Trunc = 0;
      double MaxEigen = 0;
      int MinSweeps = 1;
      int MaxSweeps = 100;
      double RealTime = 0;
      double Beta = 0;
      int NumTimesteps = 1;
      int SaveTimesteps = 0;
      std::string OutputPrefix = "";
      bool Normalize = false;
      int MaxIter = 20;
      double ConvergenceFactor = 10;
      double MixFactor = 0.001;
      double SolverTol = -1;
      bool UseK1;
      bool UseFinal = false;
      int NumPageFiles = 1;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(), 
          "operator to use for the Hamiltonian (wavefunction attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::string>(), "initial wavefunction (required)")
	 ("two-site,2", "modify 2 neighboring sites at once")
	 ("min-states,", prog_opt::value<int>(&MinStates), 
          ("Minimum number of states to keep [default " 
           +boost::lexical_cast<std::string>(MinStates)+"]").c_str())
	 ("max-states,m", prog_opt::value<int>(&MaxStates), 
          ("Maximum number of states to keep [default "
           +boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("trunc,r", prog_opt::value(&Trunc), 
          ("Desired truncation error [default "
           +boost::lexical_cast<std::string>(Trunc)+"]").c_str())
         ("eigen-cutoff,d", prog_opt::value(&MaxEigen), 
          ("Cutoff threshold for density matrix eigenvalues (alternative to truncation error) [default "
           +boost::lexical_cast<std::string>(MaxEigen)+"]").c_str())
         ("solver-tol", prog_opt::value(&SolverTol),
          "tolerance for exponential solver [default ???]")
         ("min-sweeps", prog_opt::value(&MinSweeps),
          "Minimum number of sweeps per timestep")
         ("max-sweeps", prog_opt::value(&MaxSweeps),
          "Maximum number of sweeps per timestep")
         ("convergence-factor,g", prog_opt::value(&ConvergenceFactor),
          ("Convergence criteria is |Psi-PsiPrev|^2 <= Trunc * ConvergenceFactor. "
           "default ["+boost::lexical_cast<std::string>(ConvergenceFactor)+"]").c_str())
         ("use-final", prog_opt::bool_switch(&UseFinal),
          "Do a final sweep with zero mixing immediately before advancing the simulation time")
         ("max-iter,i", prog_opt::value(&MaxIter),
          ("Minimum number of Lanczos iterations in the time evolution [default "
           +boost::lexical_cast<std::string>(MaxIter)+"]").c_str())
         ("mix-factor,f", prog_opt::value(&MixFactor), 
          ("Mixing coefficient for the density matrix [default "
           +boost::lexical_cast<std::string>(MixFactor)+"]").c_str())
	 ("timestep,t", prog_opt::value<std::complex<double> >(&Timestep), 
          "Timestep (can be complex) [default 1]")
         ("num-timesteps,n", prog_opt::value<int>(&NumTimesteps), 
          "Total number of timesteps to calculate [default 1]")
         ("save-timesteps,s", prog_opt::value<int>(&SaveTimesteps), 
          "Save the wavefunction after every s timesteps [default 0]")
         ("Time,T", prog_opt::value<double>(&RealTime), 
          "Absolute real time of the input wavefunction (wavefunction attribute \"Time\")")
	 ("Beta,B", prog_opt::value<double>(&Beta),
	  "Absolute imaginary time (inverse temperature) of the input wavefunction "
	  "(wavefunction attribute \"Beta\")")
         ("use-k1", "Construct the matrix elements <psi(t+Delta t)|H|psi(t)>")
         ("out,o", prog_opt::value<std::string>(&OutputPrefix),
          "Filename prefix for saved wavefunctions [defaults to the initial wavefunction name]")
         ("stripe-files", prog_opt::value(&NumPageFiles), 
          "Stripe the temporary disk space over this many files [default 1]")
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
         std::cerr << "usage: mp-evolve-krylov-simple [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "By default, the input wavefunction is evolved by the specified "
            "number of timesteps.\n";
         std::cerr << "If the --save-timesteps option is used, the initial wavefunction "
            "is not modified, but instead the wavefunction is saved every N steps.\n";
         return 1;
      }

      std::cout << "Starting Krylov...\n";
      std::string InputWavefunction = vm["wavefunction"].as<std::string>();
      std::cout << "Input wavefunction: " << InputWavefunction << std::endl;

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
      else
      {
         std::cout << "Evolving the input wavefunction.\n";
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
         int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
         std::string OutputDir = getenv_or_default("MP_BINPATH", std::string(""));
         if (!OutputDir.empty() && OutputDir[OutputDir.size()-1] != '/')
            OutputDir += '/';
         pheap::Initialize(OutputDir+OutputPrefix+".bin", NumPageFiles, PageSize, CacheSize);
         P = pheap::ImportHeap(InputWavefunction);
      }
      CenterWavefunction Psi = *P;
      P = pvalue_ptr<MPWavefunction>();

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
      SplitOperator Hamiltonian = ParseOperator(HamString);

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

      std::cout << "Timestep: " << Timestep << std::endl;
      // Flip the timestep to the usual form
      Timestep = complex(-Timestep.imag(), Timestep.real());

      // Clear the existing Psi attributes
      Psi.Attributes() = AttributeList();
      // And set the attributes that are relevant
      Psi.Attributes()["Hamiltonian"] = HamString;
      Psi.Attributes()["Beta"] = Beta;
      Psi.Attributes()["Time"] = RealTime;

      UseK1 = vm.count("use-k1");

      // Now we can construct the actual KrylovLoop object
      KrylovLoop dmrg(Psi, Hamiltonian, Timestep, UseK1);

      TwoSite = vm.count("two-site");
      if (TwoSite)
	 std::cout << "Optimizing two sites at a time\n";
      dmrg.SetTwoSite(TwoSite);
      std::cout << "minimum states kept: " << MinStates << '\n';
      dmrg.SetMinStates(MinStates);
      std::cout << "maximum states kept: " << MaxStates << '\n';
      dmrg.SetMaxStates(MaxStates);
      std::cout << "required truncation per site: " << Trunc << '\n';
      dmrg.SetMinTrunc(Trunc);
      std::cout << "cutoff threshold for density matrix eigenvalues: " << MaxEigen << '\n';
      dmrg.SetMaxEigen(MaxEigen);
      std::cout << "Mix factor: " << MixFactor << '\n';
      dmrg.SetMixFactor(MixFactor);
      if (SolverTol == -1)
      {
         SolverTol = MaxEigen;
      }
      std::cout << "tolerance for solver: " << SolverTol << '\n';
      dmrg.SetSolverTol(SolverTol);
      std::cout << "maximum number of Lanczos iterations: " << MaxIter << '\n';
      dmrg.SetMaxIter(MaxIter);
      std::cout << "Minimum number of sweeps per timestep: " << MinSweeps << '\n';
      dmrg.SetMinSweeps(MinSweeps);
      std::cout << "Maximum number of sweeps per timestep: " << MaxSweeps << '\n';
      dmrg.SetMaxSweeps(MaxSweeps);
      std::cout << "Convergence factor: " << ConvergenceFactor << '\n';
      dmrg.SetConvergenceFactor(ConvergenceFactor);
      std::cout << "Using final zero mixing sweep: " << (UseFinal ? "yes" : "no") << '\n';
      dmrg.SetUseFinalSweep(UseFinal);
      if (UseK1)
         std::cout << "Constructing matrix elements <psi(t+Delta t)|H|psi(t)>\n";

      dmrg.SetupLogFiles(OutputPrefix);

      for (int CurrentTimestep = 0; CurrentTimestep < NumTimesteps; ++CurrentTimestep)
      {
         dmrg.DoSweepUntilConverged();
         std::cout << "Timestep completed.  Current time is "
                   << (complex(RealTime,Beta) + complex(0,-CurrentTimestep-1)*Timestep)
                   << ", wavefunction difference is " << dmrg.Difference() << std::endl;

         if (SaveTimesteps > 0 && (CurrentTimestep+1) % SaveTimesteps == 0)
         {
            std::ostringstream filetime;
            filetime.precision(5);
            pvalue_ptr<MPWavefunction> OutPsi(new MPWavefunction(dmrg.Wavefunction().AsLinearWavefunction()));
            filetime << OutPsi->Attributes()["Time"].as<double>(); //(CurrentTimestep+1)*abs(Timestep);
            std::string FName = OutputPrefix + ".t" + filetime.str();
            std::cout << "Saving wavefunction as " << FName << '\n';
            if (Normalize)
               OutPsi.mutate()->normalize();
            pheap::ExportHeap(FName, OutPsi);
         }

         if (CurrentTimestep < NumTimesteps-1)
            dmrg.NextTimestep();
      }

      if (SaveTimesteps == 0)
      {
         P = pvalue_ptr<MPWavefunction>(new MPWavefunction(dmrg.Wavefunction().AsLinearWavefunction()));
         if (Normalize)
            P.mutate()->normalize();

         //pheap::ExportHeap("testfile", P);

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
