// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/haldane-su2.cpp
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

// SU(2) Haldane model.
//
// Example for a width-3 lattice (site numbers in brackets are periodic repeats
// in the vertical direction (i.e. top-left (5) is the same site as the
// bottom-left 5).
// Sites 6,7,8,9,10,11 are the second unit cell, e.g. 6 is (1)[0].
//
// (4)-(11)    12
//  /    \    /
//(5)     6--13
//  \    /    \
//   0--7      14
//  /    \    /
// 1      8--15
//  \    /    \
//   2--9      16
//  /    \    /
// 3      10-17
//  \    /    \
//   4--11    (12)
//  /    \    /
// 5     (6)-(13)
//  \    /
//  (0)-(7)
//
//     Y
//    /
// X--
//    \
//     Z
// Nearest-neighbour interactions are:
// (vertical)    (0)[i] (0)[(i+1)%(2w)]  (for i = 0..2w), 2w terms per unit cell
// (horizontal)  (0)[i] (1)[(i+1)%(2w)]  (for i = 0,2,4,...,2(w-1)), w terms per unit cell

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/fermion-su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

using math_const::pi;
namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int w = 4;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
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
      OpDescriptions.set_description("SU(2) Haldane model");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_t1"        , "Nearest-neighbour hopping")
         ("H_t2cw"      , "Next-nearest-neighbour clockwise hopping")
         ("H_t2acw"     , "Next-nearest-neighbour anticlockwise hopping")
         ;
      OpDescriptions.add_functions()
         ("H_t2"        , "Next-nearest-neighbour complex hopping with phase phi")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = FermionSU2();
      int u = w*2;
      UnitCell Cell = repeat(Site, u);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator I(Cell, "I"), CH(Cell, "CH"), C(Cell, "C"), Pdouble(Cell, "Pdouble"), Hu(Cell, "Hu"), N(Cell, "N");

      UnitCellMPO H_t1, H_t2cw, H_t2acw;

      for (int i = 0; i < u; i += 2)
      {
         H_t1 += dot(CH(0)[i], C(0)[(i+1)%u]) + dot(C(0)[i], CH(0)[(i+1)%u])
               + dot(CH(0)[i], C(1)[(i+1)%u]) + dot(C(0)[i], CH(1)[(i+1)%u])
               + dot(CH(0)[(i+1)%u], C(0)[(i+2)%u]) + dot(C(0)[(i+1)%u], CH(0)[(i+2)%u]);
         H_t2cw += dot(CH(0)[i], C(0)[(i+2)%u])
                 + dot(CH(0)[(i+1)%u], C(1)[(i+3)%u])
                 + dot(CH(0)[(i+2)%u], C(1)[(i+2)%u])
                 + dot(CH(1)[(i+3)%u], C(1)[(i+1)%u])
                 + dot(CH(1)[(i+2)%u], C(0)[i])
                 + dot(CH(1)[(i+1)%u], C(0)[(i+1)%u]);
         H_t2acw += dot(C(0)[i], CH(0)[(i+2)%u])
                  + dot(C(0)[(i+1)%u], CH(1)[(i+3)%u])
                  + dot(C(0)[(i+2)%u], CH(1)[(i+2)%u])
                  + dot(C(1)[(i+3)%u], CH(1)[(i+1)%u])
                  + dot(C(1)[(i+2)%u], CH(0)[i])
                  + dot(C(1)[(i+1)%u], CH(0)[(i+1)%u]);
      }

      Lattice["H_t1"] = sum_unit(H_t1);
      Lattice["H_t2cw"] = sum_unit(H_t2cw);
      Lattice["H_t2acw"] = sum_unit(H_t2acw);
      Lattice.func("H_t2")(arg("phi")) = "exp(i*phi)*H_t2cw + exp(-i*phi)*H_t2acw";

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
