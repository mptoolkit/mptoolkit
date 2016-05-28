// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-evolve-mixed.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "mp-algorithms/aggregator.h"
#include "mp-algorithms/lanczos-exponential.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/dmrg.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/proccontrol.h"
#include <boost/program_options.hpp>
#include <iostream>

namespace prog_opt = boost::program_options;

bool ShowStates = false; // for verbose output

double MinTrunc = std::numeric_limits<double>::epsilon()*8;
int MaxStates = 100000;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(SimpleOperator const& Ham_,
		      MPStateComponent const& Left_,
		      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Ham, Left, Psi, herm(Right));
   }

   SimpleOperator Ham;
   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(SimpleOperator const& Ham_,
				       MPStateComponent const& Left_,
				       MPStateComponent const& Right_)
   : Ham(Ham_), Left(Left_), Right(Right_)
{
}

int main(int argc, char** argv)
{
   try
   {
      int Location = -1;
      std::vector<std::string> Expectations;
      std::complex<double> Timestep = 1;
      int NumTimesteps = 1;
      int NumIterations = 10;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Hamiltonian,H", prog_opt::value<std::string>(), 
          "operator to use for the Hamiltonian (right Lanczos vector attribute \"Hamiltonian\")")
         ("wavefunction,w", prog_opt::value<std::vector<std::string> >(),
          "input wavefunction to generate the effective basis (one or more)")
	 ("max-states,m", prog_opt::value<int>(&MaxStates), 
          ("Maximum number of states to keep in the effective basis [default "
           + boost::lexical_cast<std::string>(MaxStates)+"]").c_str())
         ("min-trunc", prog_opt::value<double>(&MinTrunc), 
          ("Minimum desired truncation error per site of the effective basis [default "
           + boost::lexical_cast<std::string>(MinTrunc) + "]").c_str())
         ("bond,b", prog_opt::value<int>(&Location),
          "Generate the basis at this bond, valid is 1 .. L-1 [default L/2]")
	 ("timestep,t", prog_opt::value<std::complex<double> >(&Timestep), 
          "Timestep (can be complex) [default 1]")
         ("num-timesteps,n", prog_opt::value<int>(&NumTimesteps), 
          "Total number of timesteps to calculate [default 1]")
         ("expect,e", prog_opt::value<std::vector<std::string> >(&Expectations),
          "Expectation values to calculate at each timestep (zero or more)")
         ("krylov,k", prog_opt::value<int>(&NumIterations),
          "Dimension of the Krylov subspace")
         ("verbose,v", prog_opt::value<int>()->default_value(0),
          "increase verbosity")
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
         std::cerr << "usage: mp-evolve-mixed [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      int const Verbosity = vm["verbose"].as<int>();
      ShowStates = Verbosity >= 1;

      std::vector<std::string> InputWavefunctions;
      if (vm.count("wavefunction") != 0)
         InputWavefunctions = vm["wavefunction"].as<std::vector<std::string> >();

      // Load the wavefunctions
      std::string TempFile = getenv_or_default("MP_BINPATH", std::string());
      int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      if (Verbosity >= 2)
      {
         std::cerr << "Using page size " << PageSize 
                   << ", cache size " << CacheSize << '\n';
      }
      int TempFileDesc = ProcControl::CreateUniqueFile(TempFile);
      CHECK(TempFileDesc != -1);  // if this happens, then we failed to create a temporary file
      pheap::Initialize(TempFile, 1, PageSize, CacheSize);

      std::vector<MPWavefunction> Psi;
      for (unsigned i = 0; i < InputWavefunctions.size(); ++i)
      {
         if (Verbosity >= 1)
            std::cerr << "Loading wavefunction: " << InputWavefunctions[i] << '\n';
         pvalue_ptr<MPWavefunction> P = pheap::ImportHeap(InputWavefunctions[i]);
         Psi.push_back(*P);
      }

      // Set up the Hamiltonian
      std::string HamString;
      if (vm.count("Hamiltonian") == 1)
      {
         HamString = vm["Hamiltonian"].as<std::string>();
      }
      else
      {
         if (Psi[0].Attributes().count("Hamiltonian") == 0)
         {
            std::cerr << "fatal: No Hamiltonian specified.\n"
               "Specify a Hamiltonian either with --Hamiltonian parameter, "
               "or as an attribute of the right Lanczos vector.\n";
            return 1;
         }
         HamString = Psi[0].Attributes()["Hamiltonian"].as<std::string>();
      }
      if (Verbosity >= 1)
         std::cerr << "Hamiltonian: " << HamString << std::endl;

      OperatorList Lattice;
      std::vector<MPOperator> Operators;
      Operators.push_back(MPOperator());
      std::tie(Lattice, Operators[0]) = ParseLatticeAndOperator(HamString);

      // Set up the location bond
      if (Location == -1)
      {
         Location = Psi[0].size() / 2;
      }

      // Load the remaining operators
      std::vector<std::string> OperatorNames = vm.count("expect") ? 
         vm["expect"].as<std::vector<std::string> >() : std::vector<std::string>();
      for (unsigned i = 0; i < OperatorNames.size(); ++i)
      {
         Operators.push_back(Lattice[OperatorNames[i]]);
      }

      // Set up the std::cout flags
      std::cout.precision(16);
      std::cout.setf(std::ios::showpoint);
      std::cout.setf(std::ios::fixed);

      Aggregator Ag(Psi, Operators, ShowStates, MaxStates, MinTrunc, Location);

      std::complex<double> Theta(Timestep.imag(), Timestep.real());
      TRACE(Theta);

      MatrixOperator State = Ag.Center(0);
      std::cout << 0.0 << "  ";
      for (unsigned j = 0; j < OperatorNames.size(); ++j)
      {
         std::cout << Ag.Expectation(State, j+1) << "   ";
      }
      std::cout << std::endl;
      for (int i = 0; i < NumTimesteps; ++i)
      {
         State = LanczosExponential(State, 
                                    SuperblockMultiply(conj(Ag.OpCenter(0)),
                                                       Ag.OpLeft(0),
                                                       Ag.OpRight(0)), 
                                    NumIterations, 
                                    Theta);

         std::cout << (Timestep.real() * i) << "  ";
         for (unsigned j = 0; j < OperatorNames.size(); ++j)
         {
            std::cout << Ag.Expectation(State, j+1) << "   ";
         }
         std::cout << std::endl;
      }

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
