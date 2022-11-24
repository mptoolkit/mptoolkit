// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-finegrain.cpp
//
// Copyright (C) 2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "lattice/infinitelattice.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

template <typename T>
std::vector<T>
repeat(std::vector<T> v, int N)
{
   std::vector<T> Result;
   Result.reserve(v.size() * N);
   for (int i = 0; i < N; ++i)
   {
      for (auto const& x : v)
      {
         Result.push_back(x);
      }
   }
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OutputFile;
      std::string InputFile;
      bool Force = false;
      int Finegrain = 2;
      std::string LatticeFile;
      int MinStates = 1;
      int States = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("lattice,l", prog_opt::value(&LatticeFile), "lattice file for the fine-grained lattice [required]")
         ("force,f", prog_opt::bool_switch(&Force),
          "allow overwriting the output file, if it already exists")
	 ("finegrain", prog_opt::value(&Finegrain),
	  FormatDefault("coarse-grain i-to-N", Finegrain).c_str())
         ("min-states", prog_opt::value<int>(&MinStates),
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("states,m", prog_opt::value(&States),
          FormatDefault("number of states", States).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") < 1 || vm.count("lattice") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -l <lattice> <input_psi> <output_psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = States;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFile == OutputFile)
         PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFile);
      }

      // load the wavefunction
      InfiniteWavefunctionLeft PsiL = PsiPtr->get<InfiniteWavefunctionLeft>();

      LinearWavefunction Psi;
      RealDiagonalOperator Lambda;
      std::tie(Psi, Lambda) = get_left_canonical(PsiL);
      QuantumNumbers::QuantumNumber QShift = PsiPtr->get<InfiniteWavefunctionLeft>().qshift();

      // load the lattice and get the list of basis
      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);
      std::vector<BasisList> FullBasis = Basis1FromSiteList(*Lattice->GetUnitCell().GetSiteList());
      int UnitCellSize = Psi.size() * Finegrain;
      if (UnitCellSize % FullBasis.size() != 0)
      {
         std::cerr << "Wavefunction size is not compatible with the fine-grained lattice.\n";
         return 1;
      }
      FullBasis = repeat(FullBasis, UnitCellSize / FullBasis.size());

      // Do the fine-grain operation
      MatrixOperator M = Lambda;
      std::tie(M, Psi) = fine_grain(Psi, M, FullBasis, Finegrain, SInfo, Verbose);

      // convert back to an InfiniteWavefunctionLeft. It is better to do this by left orthogonalizing
      // again, since we're in completely the wrong basis now
      M = left_orthogonalize(M, Psi, Verbose-1);
      // And shift to the basis where Lambda is diagonal
      MatrixOperator U, Vh;
      SingularValueDecomposition(M, U, Lambda, Vh);
      Psi.set_back(prod(Psi.get_back(), U*Vh));
      M = Lambda;

      M = herm(M) * M; // convert to a density matrix
      PsiL = InfiniteWavefunctionLeft::Construct(Psi, QShift, M, PsiL.log_amplitude()/Finegrain, Verbose-1);
      Psi = LinearWavefunction(); // destroy it, we don't need it anymore

      PsiPtr.mutate()->Wavefunction() = PsiL;
      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
