// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spinlessfermionladder-u1u1.cpp
//
// Copyright (C) 2012-2023 Ian McCulloch <ian@qusim.net>
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

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.description("Spinless fermion ladder U(1)xU(1) symmetry");
      OpDescriptions.author("IP McCulloch", "ian@qusim.net");
      OpDescriptions.add_operators()
         ("H_t1"  , "nearest neighbor hopping in the first leg")
         ("H_t2"  , "nearest neighbor hopping in the second leg")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      LatticeSite Site1 = SpinlessFermionU1("N1");
      LatticeSite Site2 = SpinlessFermionU1("N2");
      UnitCell Cell(QuantumNumbers::SymmetryList("N1:U(1),N2:U(1)"), Site1, Site2);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N");

      Lattice["H_t1"]  = -sum_unit(dot(CH(0)[0], C(1)[0]) - dot(C(0)[0], CH(1)[0]));
      Lattice["H_t2"]  = -sum_unit(dot(CH(0)[1], C(1)[1]) - dot(C(0)[1], CH(1)[1]));

      // Information about the lattice
      Lattice.set_description("U(1)x(1) Spinless Fermion Fermi-Hubbard 2-leg ladder");
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice
      pheap::ExportObject(FileName, Lattice);
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
