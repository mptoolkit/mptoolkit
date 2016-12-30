// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-rename-symmetry.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "interface/inittemp.h"
#include "interface/operator-parser.h"
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

LinearWavefunction RenameSymmetry(LinearWavefunction const& Psi, SymmetryList const& NewSL)
{
   LinearWavefunction Result(NewSL);
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(RenameSymmetry(*I, NewSL));
   }
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      if (argc != 4)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-rename-symmetry <old> <new> <psi>\n";
         return 1;
      }

      std::string OldQ = argv[1];
      std::string NewQ = argv[2];
      std::string FName = argv[3];
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

      int const s = Psi->GetSymmetryList().WhichSymmetry(OldQ);
      if (s == -1)
      {
         std::cerr << "mp-shift-basis: warning: \"" + OldQ
            + "\" does not name a quantum number.\n";
         pheap::ShutdownPersistent(Psi);
         return 0; // this isn't really an error
      }

      // construct the new symmetry list
      SymmetryList NewSL = ChangeSymmetryName(Psi->GetSymmetryList(), OldQ, NewQ);

      Psi = pvalue_ptr<MPWavefunction>(new MPWavefunction(RenameSymmetry(*Psi, NewSL)));
      pheap::ShutdownPersistent(Psi);

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
