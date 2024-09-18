// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spincylinder-su2.cpp
//
// Copyright (C) 2015-2023 Ian McCulloch <ian@qusim.net>
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
#include "models/spin-su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int x = 0;
      int y = 4;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the spin", Spin).c_str())
         (",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str())
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector", y).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.description("SU(2) spin cylinder");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1x"  , "nearest neighbor spin exchange in the x direction")
         ("H_J1y"  , "nearest neighbor spin exchange in the y direction")
         ("H_J1"   , "nearest neighbor spin exchange H_J1x + H_J2x")
         ("H_J2d"  , "next-nearest neighbor spin exchange in the diagonal (0,0) - (1,1) direction")
         ("H_J2a"  , "next-nearest neighbor spin exchange in the antidiagonal (0,1) - (1,0) direction")
         ("H_J2"   , "next-nearest neighbor spin exchange H_J1a + H_J2d")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Spin cylinder.  Wrapping vector is (x,y).  y is the width, x is the offset.  (0,y) is the\n"
                   << "YC configuration with width and unit cell size y.  (x,y) is with y cylinder\n"
                   << "with offset x.  The unit cell size is x.\n";
         std::cout << OpDescriptions << '\n';
         return 1;
      }

      int CellSize = x == 0 ? y : x;

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator S(Cell, "S");

      UnitCellMPO J1x, J1y, J2d, J2a;
      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            J1x += inner(S(0)[i], S(1)[i]);
            J1y += inner(S(0)[i], S(0)[(i+1)%y]);
            J2d += inner(S(0)[i], S(1)[(i+1)%y]);
            J2a += inner(S(0)[(i+1)%y], S(1)[i]);
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            J1x += inner(S(0)[i], S(0)[i+1]);
            J2d += inner(S(0)[i], S(1)[i+1]);
         }
         J1x += inner(S(0)[x-1], S(y+1)[0]);
         J2d += inner(S(0)[x-1], S(y+2)[0]);
         for (int i = 0; i < x; ++i)
         {
            J1y += inner(S(0)[i], S(1)[i]);
            J2a += inner(S(0)[(i+1)%x], S(1)[i]);
         }
      }

      Lattice["H_J1x"] = sum_unit(J1x);
      Lattice["H_J1y"] = sum_unit(J1y);
      Lattice["H_J1"] = sum_unit(J1x+J1y);
      Lattice["H_J2d"] = sum_unit(J2d);
      Lattice["H_J2a"] = sum_unit(J2a);
      Lattice["H_J2"] = sum_unit(J2d+J2a);

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
