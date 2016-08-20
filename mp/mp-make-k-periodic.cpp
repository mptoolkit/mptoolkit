// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-make-k-periodic.cpp
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
#include "matrixproduct/operatoratsite.h"
#include "mp/copyright.h"
#include "common/math_const.h"
#include "pheap/pheap.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <cmath>
#include "common/environment.h"

namespace prog_opt = boost::program_options;
using std::sin;

int main(int argc, char** argv)
{
   try
   {
      std::string LatticeFile;
      std::string OperatorName;
      std::string OutputName;
      int UnitCellSize = 1;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("lattice,l", prog_opt::value(&LatticeFile),
          "lattice file [required]")
         ("input-operator,i", prog_opt::value(&OperatorName),
          "name of the input operator [required]")
         ("output-operator,o", prog_opt::value(&OutputName),
          "name of the output operator [defaults to <input-name>_k]")
         ("unit-size,u", prog_opt::value(&UnitCellSize),
          "size of the unit cell [default 1]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || !vm.count("lattice") || !vm.count("input-operator"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-make-k [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (!vm.count("output-operator"))
         OutputName = OperatorName + "_k";

      // load the lattice file
      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<OperatorList> System = pheap::OpenPersistent(LatticeFile, CacheSize);

      // get the lattice size, accounting for the unit cell
      int LatticeSize = System->size();
      CHECK(LatticeSize % UnitCellSize == 0)("The total lattice size is not a multiple of the unit cell!");
      LatticeSize /= UnitCellSize;

      OperatorAtSite<OperatorList const, int> InputOp(*System, OperatorName);
      OperatorAtSite<OperatorList, int> OutputOp(*System.mutate(), OutputName);

      for (int k = -LatticeSize+1; k <= LatticeSize; ++k)
      {
         std::cout << "k=" << k << '\n';
         MPOperator OpK;
         for (int i = 0; i < LatticeSize; ++i)
         {
            OpK += std::sqrt(2./double(LatticeSize)) *
               std::complex<double>(cos(math_const::pi * k * i / LatticeSize),
                                    sin(math_const::pi * k * i / LatticeSize))
               * InputOp(i+1);
         }
         OutputOp(k) = OpK;
      }

      // save the lattice
      pheap::ShutdownPersistent(System);

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
