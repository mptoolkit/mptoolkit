// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spinlessfermion-u1.cpp
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
      OpDescriptions.description("Spinless fermion U(1) symmetry");
      OpDescriptions.author("J Pillay", "pillayjason@hotmail.com");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_t"   , "nearest neighbor hopping")
         ("H_tc"  , "complex nearest neighbor hopping")
         ("H_U"   , "nearest-, next-nearest- and next-next-nearest neighbour Coulomb interaction")
         ;

      // Descriptions for the operators
      OpDescriptions.add_functions()
         ("H_V", "3-site periodic potential, parametized by angle 'delta'")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      LatticeSite Site = SpinlessFermionU1();
      UnitCell Cell(Site);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N");

      Lattice["H_t"]  = sum_unit(dot(CH(0), C(1)) - dot(C(0), CH(1)));
      Lattice.func("H_t3")(arg("theta")) = "sum_unit(sites=3, dot(CH(0), C(1)) - dot(C(0), CH(1))"
                                                          " + dot(CH(1), C(2)) - dot(C(1), CH(2))"
                                                          " + exp(i*theta)*dot(CH(2), C(3)) - exp(-i*theta)*dot(C(2), CH(3)))";
      Lattice["H_tc"] = sum_unit(std::complex<double>(0,1)*dot(CH(0), C(1)) - dot(C(0), CH(1)));
      Lattice["H_U"]  = sum_unit(dot(N(0), N(1))) + sum_unit(dot(N(0), N(2))) + sum_unit(dot(N(0), N(3)));

      Lattice.func("H_V")(arg("delta")) = "sum_unit(sites=3, cos(delta)*N(0) + cos(2*pi/3 + delta)*N(1) + cos(4*pi/3 + delta)*N(2))";

      // Information about the lattice
      Lattice.set_description("U(1) Spinless Fermion Fermi-Hubbard model");
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
