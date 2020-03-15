// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/toric-code-2.cpp
//
// Copyright (C) 2020 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// Toric code model with a two-site unit cell.
//
// The lattice for a cylinder with w=2 is shown below (bracketed numbers are
// periodically repeated sites).
//
//             +-(-1)-+-(3)-+-(7)-+
//             |     |     |
//             0     4     8
//             |     |     |
//       +--1--+--5--+--9--+
//       |     |     |     
//       2     6     10
//       |     |     |
// +--3--+--7--+--11-+
//
// Alternatively, rotating by 45 degrees so the unit cell is vertical.
//
// \    /\
//(-2)(3) 8
//   \/    \/
//   /\    /\
//(-1) 4  9
// /    \/
// \    /\
//  0  5  10
//   \/    \/
//   /\    /\
//  1  6  11
// /    \/
// \    /\
//  2  7 (12)
//   \/    \/
//   /\    /\
//  3 (8)(13)
// /    \/

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 4;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the spin", Spin).c_str())
         ("width,w", prog_opt::value(&w), FormatDefault("width of the cylinder", w).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.description("Toric code with a two-site unit cell");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("A"          , "star operator")
         ("B"          , "plaquette operator")
         ;
      OpDescriptions.add_operators()
         ("H_x"        , "magnetic field in x direction")
         ("H_y"        , "magnetic field in y direction")
         ("H_z"        , "magnetic field in z direction")
         ("H_star"     , "sum of star operators")
         ("H_plaq"     , "sum of plaquette operators")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cout << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = SpinSite(Spin);
      UnitCell Cell(repeat(Site, 2));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator I(Cell, "I"), X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");
      UnitCellOperator A(Cell, "A"), B(Cell, "B");
      
      // Magnetic fields.
      Lattice["H_x"] = sum_unit(X(0)[0] + X(0)[1]);
      Lattice["H_y"] = sum_unit(Y(0)[0] + Y(0)[1]);
      Lattice["H_z"] = sum_unit(Z(0)[0] + Z(0)[1]);

      // Star and plaquette operators.
      A = X(0)[0] * X(0)[1] * X(w)[1] * X(w+1)[0];
      B = Z(0)[1] * Z(1)[0] * Z(w+1)[0] * Z(w+1)[1];

      // Sum of the star and plaquette operators.
      Lattice["H_star"] = sum_unit(A);
      Lattice["H_plaq"] = sum_unit(B);

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
