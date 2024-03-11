// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-scale-basis.cpp
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

LinearWavefunction ScaleBasisU1(LinearWavefunction const& Psi,
                                     std::string const& Name,
                                     double Factor)
{
   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(ScaleBasisU1(*I, Name, Factor));
   }
   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-scale-basis <quantum-number-name> <factor> <psi>\n";
      std::cerr << "this is applicable only to U(1) quantum numbers.\n";
      return 1;
   }

   std::string QName = argv[1];
   double Factor = boost::lexical_cast<double>(argv[2]);
   std::string FName = argv[3];
   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

   // verify that QName is actually a symmetry
   int const s = Psi->GetSymmetryList().WhichSymmetry(QName);
   if (s == -1)
   {
      std::cerr << "mp-scale-basis: warning: \"" + QName
         + "\" does not name a quantum number.\n";
      pheap::ShutdownPersistent(Psi);
      return 0; // this isn't really an error
   }
   // else
   if (Psi->GetSymmetryList().SymmetryType(s) != "U(1)")
   {
      std::cerr << "mp-scale-basis: error: \"" + QName +
         "\" is not a U(1) symmetry, no action taken.\n";
      pheap::ShutdownPersistent(Psi);
      return 2;
   }

   Psi = pvalue_ptr<MPWavefunction>(new MPWavefunction(ScaleBasisU1(*Psi, QName, Factor)));
   pheap::ShutdownPersistent(Psi);
}
