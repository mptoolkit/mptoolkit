// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinladder-su2.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int Legs = 2;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the spin", Spin).c_str())
         ("legs,l", prog_opt::value(&Legs), FormatDefault("number of legs", Legs).c_str())
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
      OpDescriptions.description("SU(2) spin ladder");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1x"  , "nearest neighbor spin exchange in the x (infinite) direction")
         ("H_J1y"  , "nearest neighbor spin exchange in the y (rung) direction")
         ("H_J1yp" , "nearest neighbor spin exchange in the y (rung) direction, including periodic term")
         ("H_J1"   , "nearest neighbor spin exchange H_J1x + H_J2x")
         ("H_J1p"  , "nearest neighbor spin exchange H_J1xp + H_J2x (periodic in y direction)")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Square lattice spin ladder.  X direction is the long (infinite) direction, \n"
                   << "Y direction is the short direction.\n";
         std::cerr << OpDescriptions << '\n';
            ;
         return 1;
      }

      int CellSize = Legs;

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator S(Cell, "S");

      UnitCellMPO J1x, J1y;
      for (int i = 0; i < Legs-1; ++i)
      {
         J1x += inner(S(0)[i], S(1)[i]);
         J1y += inner(S(0)[i], S(0)[i+1]);
      }
      J1x += inner(S(0)[Legs-1], S(1)[Legs-1]);

      Lattice["H_J1x"] = sum_unit(J1x);
      Lattice["H_J1y"] = sum_unit(J1y);
      Lattice["H_J1yp"] = sum_unit(J1y + inner(S(0)[0], S(0)[Legs-1]));
      Lattice["H_J1"] = sum_unit(J1x+J1y);
      Lattice["H_J1p"] = sum_unit(J1x+J1y + inner(S(0)[0], S(0)[Legs-1]));

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
