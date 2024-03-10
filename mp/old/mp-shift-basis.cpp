// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-shift-basis.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
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

LinearWavefunction ShiftLocalBasis(LinearWavefunction const& Psi, QuantumNumber LocalShift)
{
   // the running quantum number that we shift the matrix basis.  initially zero
   QuantumNumber MatrixShift(Psi.GetSymmetryList());

   LinearWavefunction Result(Psi.GetSymmetryList());
   LinearWavefunction::const_iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;
      Result.push_front(ShiftLocalBasis(*I, LocalShift, MatrixShift));
      MatrixShift = transform_targets(MatrixShift, LocalShift)[0]; // update running shift
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
         std::cerr << "usage: mp-shift-basis <quantum-number-name> <shift> <psi>\n";
         return 1;
      }

      std::string QName = argv[1];
      std::string Shift = argv[2];
      std::string FName = argv[3];
      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize());

      // get the symmetry type
      int const s = Psi->GetSymmetryList().WhichSymmetry(QName);
      if (s == -1)
      {
         std::cerr << "mp-shift-basis: warning: \"" + QName
            + "\" does not name a quantum number.\n";
         pheap::ShutdownPersistent(Psi);
         return 0; // this isn't really an error
      }

      // construct the shift quantum number
      // firstly, a symmetry list consisting only of the relevant symmetry
      std::string SL = QName + ':' + Psi->GetSymmetryList().SymmetryType(s);
      QuantumNumber ShiftQ(SL, Shift); // construct the 1-component quantum number
      ShiftQ.Coerce(Psi->GetSymmetryList()); // map that into the full symmetry list

      Psi = pvalue_ptr<MPWavefunction>(new MPWavefunction(ShiftLocalBasis(*Psi, ShiftQ)));
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
