// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/klm-u1su2.cpp
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

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/kondo-u1su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      // Parameters of the lattice (with defaults, if applicable)
      std::string FileName;

      // Options
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

      // Descriptions of each operator
      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1)xSU(2) Kondo Lattice Model");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_t"   , "nearest neighbor hopping")
         ("H_t2"  , "nexxt-nearest neighbor hopping")
         ("H_U"   , "on-site Coulomb interaction")
         ("H_K"   , "Kondo coupling")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = KondoU1SU2();

      // The UnitCell consists of a single site
      UnitCell Cell(Site);

      // Make an infinite lattice of our unit cell
      InfiniteLattice Lattice(&Cell);

      // A short-cut to refer to an operator defined within our unit cell
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), ScSf(Cell, "ScSf"), Pdouble(Cell, "Pdouble");

      // Define operators that have support over the infinite lattice
      Lattice["H_t"]  = -sum_unit(dot(CH(0), C(1)) + dot(C(0), CH(1)));
      Lattice["H_t2"] = -sum_unit(dot(CH(0), C(2)) + dot(C(0), CH(2)));
      Lattice["H_K"]  =  sum_unit(ScSf(0));
      Lattice["H_U"] =  sum_unit(Pdouble(0));

      // Information about the lattice
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice to disk
      pheap::ExportObject(FileName, Lattice);
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
