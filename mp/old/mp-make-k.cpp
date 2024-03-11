// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-make-k.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

std::string replace_template(std::string const& s, int n)
{
   std::string Result = s;
   int pos = Result.find("%n");
   CHECK(pos >= 0)("expecting to find '%n' somewhere in the string")(s)(pos);
   Result.replace(pos, 2, boost::lexical_cast<std::string>(n));
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      std::string LatticeFile;
      std::string OperatorName;
      std::string OutputName;
      int UnitCellSize = 1;
      std::string Template = "%n";
      bool Verbose = false;

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
         ("template,t", prog_opt::value(&Template),
          "template for the coordinates to use [default '%n']")
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "increase verbosity")
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

      for (int k = 1; k <= LatticeSize; ++k)
      {
         std::cout << "k=" << k << '\n';
         std::string CoordK = replace_template(Template, k);
         MPOperator OpK;
         for (int i = 1; i <= LatticeSize; ++i)
         {
            std::string CoordI = replace_template(Template, i);
            if (Verbose)
               std::cerr << "Adding component at site " << i << " coordinates " << CoordI << '\n';
            OpK += std::sqrt(2./double(LatticeSize+1)) * sin(math_const::pi * k * i / (LatticeSize+1)) * System->Lookup(OperatorName, CoordI);
         }
         (*System.mutate())[OutputName+"("+CoordK+")"] = OpK;
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
