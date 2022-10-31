// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spincylinder-z2.cpp
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "models/spin-z2.h"
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
      OpDescriptions.description("Z2 spin cylinder");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_xx_x"  , "nearest neighbor spin coupling Sx Sx in the x direction")
         ("H_yy_x"  , "nearest neighbor spin coupling Sy Sy in the x direction")
         ("H_zz_x"  , "nearest neighbor spin coupling Sz Sz in the x direction")
         ("H_xx_y"  , "nearest neighbor spin coupling Sx Sx in the y direction")
         ("H_yy_y"  , "nearest neighbor spin coupling Sy Sy in the y direction")
         ("H_zz_y"  , "nearest neighbor spin coupling Sz Sz in the y direction")
         ("H_xx"    , "nearest neighbor spin coupling Sx Sx")
         ("H_yy"    , "nearest neighbor spin coupling Sy Sy")
         ("H_zz"    , "nearest neighbor spin coupling Sz Sz")
         ("H_J1z"   , "same as H_zz")
         ("H_J1t"   , "transverse spin exchange, H_xx + H_yy")
         ("H_J1"    , "nearest neighbor spin exchange = H_J1z + H_J1t")
         ("H_x"     , "magnetic field in the x direction")
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

      LatticeSite Site = SpinZ2(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz"), Sp(Cell, "Sp"), Sm(Cell, "Sm");

      UnitCellMPO H_xx_x, H_yy_x, H_zz_x, H_xx_y, H_yy_y, H_zz_y;
      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            H_xx_x += Sx(0)[i]*Sx(1)[i];
            H_yy_x += Sy(0)[i]*Sy(1)[i];
            H_zz_x += Sz(0)[i]*Sz(1)[i];
            H_xx_y += Sx(0)[i]*Sx(0)[(i+1)%y];
            H_yy_y += Sy(0)[i]*Sy(0)[(i+1)%y];
            H_zz_y += Sz(0)[i]*Sz(0)[(i+1)%y];
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            H_xx_x += Sx(0)[i]*Sx(0)[i+1];
            H_yy_x += Sy(0)[i]*Sy(0)[i+1];
            H_zz_x += Sz(0)[i]*Sz(0)[i+1];
         }
         H_xx_x += Sx(0)[x-1]*Sx(y+1)[0];
         H_yy_x += Sy(0)[x-1]*Sy(y+1)[0];
         H_zz_x += Sz(0)[x-1]*Sz(y+1)[0];
         for (int i = 0; i < x; ++i)
         {
            H_xx_y += Sx(0)[i]*Sx(1)[i];
            H_yy_y += Sy(0)[i]*Sy(1)[i];
            H_zz_y += Sz(0)[i]*Sz(1)[i];
         }
      }

      Lattice["H_xx_x"] = sum_unit(H_xx_x);
      Lattice["H_yy_x"] = sum_unit(H_yy_x);
      Lattice["H_zz_x"] = sum_unit(H_zz_x);
      Lattice["H_xx_y"] = sum_unit(H_xx_y);
      Lattice["H_yy_y"] = sum_unit(H_yy_y);
      Lattice["H_zz_y"] = sum_unit(H_zz_y);
      Lattice["H_xx"] = sum_unit(H_xx_x+H_xx_y);
      Lattice["H_yy"] = sum_unit(H_yy_x+H_yy_y);
      Lattice["H_zz"] = sum_unit(H_zz_x+H_zz_y);
      Lattice["H_J1z"] = Lattice["H_zz"];
      Lattice["H_J1t"] = Lattice["H_xx"] + Lattice["H_yy"];
      Lattice["H_J1"] = Lattice["H_xx"] + Lattice["H_yy"] + Lattice["H_zz"];

      // Magnetic fields.
      UnitCellMPO H_x;

      for (int i = 0; i < CellSize; ++i)
         H_x += Sx(0)[i];

      Lattice["H_x"] = sum_unit(H_x);

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
