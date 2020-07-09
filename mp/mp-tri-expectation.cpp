// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-imoments.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "wavefunction/mpwavefunction.h"
#include "wavefunction/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "common/polynomial.h"
#include "tensor/tensor_eigen.h"
#include "mp/copyright.h"
#include "common/prog_options.h"
#include "lattice/infinite-parser.h"
#include "wavefunction/momentum_operations.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "common/prog_opt_accum.h"
#include <boost/algorithm/string.hpp>
#include "common/randutil.h"

namespace prog_opt = boost::program_options;
using formatting::format_complex;

int main(int argc, char** argv)
{
   std::string FName;
   std::string OpStr;
   int N = 100;
   bool ShowAll = false;
   bool Randomize = false;

   int Verbose = 0;

   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
	 ("sites,n", prog_opt::value(&N), "maximum length to calculate (sites)")
	 ("showall", prog_opt::bool_switch(&ShowAll), "show all columns of the MPO")
	 ("randomize", prog_opt::bool_switch(&Randomize), "randomize boundary tensors")
         ("help", "show this help message")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "extra debug output (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&FName), "wavefunction")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("operator", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " <psi1> <operator>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(FName, CacheSize, true);
      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();
      int WavefuncUnitCellSize = Psi.size();

      if (!vm.count("operator"))
      {
         OpStr = PsiPtr->Attributes()["Hamiltonian"].get_or_default(std::string());
         if (OpStr.empty())
         {
            std::cerr <<  basename(argv[0]) << ": fatal: no operator specified, and wavefunction "
               "attribute Hamiltonian does not exist or is empty.\n";
            return 1;
         }
      }

      BasicTriangularMPO Op;
      InfiniteLattice Lattice;
      std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(OpStr);

      // Check that the local basis for the wavefunction and hamiltonian are compatible
      if (!is_compatible(ExtractLocalBasis(Psi), ExtractLocalBasis1(Op)))
      {
         std::cerr << "fatal: operator is defined on a different local basis to the wavefunction.\n";
	 TRACE(ExtractLocalBasis(Psi))(ExtractLocalBasis1(Op));
         return 1;
      }

      if (ExtractLocalBasis1(Op) != ExtractLocalBasis2(Op))
      {
         std::cerr << "fatal: operator has different domain and co-domain.\n";
         return 1;
      }

      StateComponent E = Initial_E(Op, Psi.Basis1());
      if (Randomize)
      {
	 srand(randutil::crypto_rand());
	 for (unsigned j = 1; j < E.size(); ++j)
	 {
	    E[j] = MakeRandomMatrixOperator(E[j].Basis1(), E[j].Basis2(), E[j].TransformsAs());
	 }
      }
      InfiniteWavefunctionLeft::const_mps_iterator PsiI = Psi.begin();
      InfiniteWavefunctionLeft::const_lambda_iterator RhoI = Psi.lambda_begin();
      ++RhoI;
      for (unsigned i = 0; i < N; ++i)
      {
	 E = contract_from_left(Op[i%Op.size()], herm(*PsiI), E, *PsiI);
	 MatrixOperator Rho = *RhoI;
	 if (ShowAll)
	 {
	    for (unsigned j = 0; j < E.size(); ++j)
	    {
	       if (is_scalar(E[j].TransformsAs()))
	       {
		  std::cout << (i+1) << " column " << (j+1) << ' ' << format_complex(inner_prod(E[j], Rho)) << '\n';
	       }
	    }
	 }
	 else
	 {
	    if (E.back().TransformsAs() == Rho.TransformsAs())
	    {
	       std::cout << (i+1) << ' ' << format_complex(inner_prod(E.back(), Rho)) << '\n';
	    }
	 }
	 ++PsiI;
	 if (PsiI == Psi.end())
	    PsiI = Psi.begin();
	 ++RhoI;
	 if (RhoI == Psi.lambda_end())
	 {
	    RhoI = Psi.lambda_begin();
	    ++RhoI;
	 }
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
