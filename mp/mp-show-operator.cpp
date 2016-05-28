// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-show-operator.cpp
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

#include "lattice/infinite-parser.h"
#include "lattice/unitcell-parser.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::vector<std::string> FiniteOperators;
      std::vector<std::string> ProductOperators;
      std::vector<std::string> TriangularOperators;
      double UnityEpsilon = DefaultClassifyUnityEpsilon;
      int Verbose = 0;
      bool Optimize = false;
      bool CoarseGrain = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("finite,f", prog_opt::value(&FiniteOperators), "parse a finite MPO")
	 ("product,p", prog_opt::value(&ProductOperators), "parse a product MPO")
	 ("triangular,t", prog_opt::value(&TriangularOperators), "parse a triangular MPO")
	 ("coarsegrain", prog_opt::bool_switch(&CoarseGrain), 
	  "for a finite operator, go through a coarse-graining then fine-graining sequence")
	 ("optimize", prog_opt::bool_switch(&Optimize), "Optimize the operator expression")
	 ("unityepsilon", prog_opt::value(&UnityEpsilon),
	  FormatDefault("Epsilon value for testing eigenvalues near unity", UnityEpsilon).c_str())
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::positional_options_description p;

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || (FiniteOperators.empty() && ProductOperators.empty() && TriangularOperators.empty())) 
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " -f|-p|-t Operator\n";
         std::cerr << desc << "\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      mp_pheap::InitializeTempPHeap();

      for (unsigned i = 0; i < FiniteOperators.size(); ++i)
      {
	 std::cout << "Finite Operator " << FiniteOperators[i] << '\n';
	 UnitCellMPO Op;
	 InfiniteLattice Lattice;
	 std::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(FiniteOperators[i]);
	 if (CoarseGrain)
	 {
	    SimpleOperator S = coarse_grain(Op.MPO()).scalar();
	    Op.MPO() = fine_grain(S, Op.MPO().LocalBasis1List(), Op.MPO().LocalBasis2List());
	 }
	 if (Optimize)
	    optimize(Op);
	 print_structure(Op.MPO(), std::cout);
	 if (Verbose > 0)
	 {
	    std::cout << Op.MPO() << '\n';
	 }
	 if (Verbose > 1)
	 {
	    SimpleRedOperator x = coarse_grain(Op.MPO());
	    std::cout << x << "\n";
	 }
      }

      for (unsigned i = 0; i < ProductOperators.size(); ++i)
      {
	 std::cout << "Product Operator " << ProductOperators[i] << '\n';
	 ProductMPO Op;
	 InfiniteLattice Lattice;
	 std::tie(Op, Lattice) = ParseProductOperatorAndLattice(ProductOperators[i]);     
	 //	 if (!NoOptimize)
	 //	    optimize(Op);
	 print_structure(Op, std::cout, UnityEpsilon);
	 if (Verbose > 0)
	 {
	    std::cout << Op << '\n';
	 }
      }
      
      for (unsigned i = 0; i < TriangularOperators.size(); ++i)
      {
	 std::cout << "Triangular Operator " << TriangularOperators[i] << '\n';
	 TriangularMPO Op;
	 InfiniteLattice Lattice;
	 std::tie(Op, Lattice) = ParseTriangularOperatorAndLattice(TriangularOperators[i]);     
	 if (Optimize)
	    optimize(Op);
	 print_structure(Op, std::cout, UnityEpsilon);
	 if (Verbose > 0)
	 {
	    std::cout << Op << '\n';
	 }
      }
      
   }
   catch(std::exception& e) 
   {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
   }
   catch(...) 
   {
      std::cerr << "Exception of unknown type!\n";
   }
   return 0;
}
