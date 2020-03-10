// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/kitaev-hex-2-z2.cpp
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

// Z2 Kitaev honeycomb model with a two-site unit cell.
//
// Horizontal bonds are Sx, +60 degree bonds are Sy, -60 degree bonds are Sz.
//
// Example for a width-3 lattice (site numbers in brackets are periodic repeats
// in the vertical direction (i.e. top-left (5) is the same site as the
// bottom-left 5).
// 
// The lattice site n on the graph below corresponds to (n/2)[n%2].
//
//(-2)-(5)    12
//  /    \    /
//(-1)    6--13
//  \    /    \
//   0--7      14
//  /    \    /
// 1      8--15
//  \    /    \
//   2--9      16
//  /    \    /
// 3      10-17
//  \    /    \
//   4--11    (18)
//  /    \    /
// 5    (12)-(19)
//  \    /
//  (6)-(13)
//
//     Y
//    /
// X--
//    \
//     Z

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-z2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

using math_const::pi;
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

      // Descriptions of each operator
      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("Z2 Kitaev honeycomb model with a two-site unit cell");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("W"          , "Plaquette operator")
         ;
      OpDescriptions.add_operators()
         ("H_x"        , "Horizontal Sx-Sx interaction")
         ("H_y"        , "+60 degrees Sy-Sy interaction")
         ("H_z"        , "-60 degrees Sz-Sz interaction")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = SpinZ2(Spin);
      UnitCell Cell = repeat(Site, 2);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator I(Cell, "I"), X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");
      UnitCellOperator W(Cell, "W");

      // Plaquette operator.
      W = Z(0)[0] * X(0)[1] * Y(1)[0] * Z(w+1)[1] * X(w+1)[0] * Y(w)[1];

      // The operators have a minus sign, following Kitaev convention
      Lattice["H_x"] = -sum_unit(inner(X(0)[0], X(w)[1]));
      Lattice["H_y"] = -sum_unit(inner(Y(0)[0], Y(0)[1]));
      Lattice["H_z"] = -sum_unit(inner(Z(0)[0], Z(-1)[1]));

      // Information about the lattice
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice to disk
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
