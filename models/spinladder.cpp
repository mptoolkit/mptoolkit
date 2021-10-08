// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinladder-su2.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "models/spin.h"
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
      OpDescriptions.description("spin ladder");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1_x"  , "nearest neighbor spin exchange in the x (infinite) direction")
         ("H_J1_y"  , "nearest neighbor spin exchange in the y (rung) direction")
         ("H_J1_yp" , "nearest neighbor spin exchange in the y (rung) direction, including periodic term")
         ("H_J1"    , "nearest neighbor spin exchange H_J1x + H_J2x")
         ("H_J1_p"  , "nearest neighbor spin exchange H_J1xp + H_J2x (periodic in y direction)")
         ("H_xx_x"  , "nearest neighbor spin coupling Sx Sx in the x direction")
         ("H_yy_x"  , "nearest neighbor spin exchange Sy Sy in the x direction")
         ("H_zz_x"  , "nearest neighbor spin exchange Sz Sz in the x direction")
         ("H_xx_y"  , "nearest neighbor spin coupling Sx Sx in the y direction")
         ("H_yy_y"  , "nearest neighbor spin exchange Sy Sy in the y direction")
         ("H_zz_y"  , "nearest neighbor spin exchange Sz Sz in the y direction")
         ("H_xx_yp" , "nearest neighbor spin coupling Sx Sx in the periodic y direction")
         ("H_yy_yp" , "nearest neighbor spin exchange Sy Sy in the periodic y direction")
         ("H_zz_yp" , "nearest neighbor spin exchange Sz Sz in the periodic y direction")
         ("H_x"     , "magnetic field in the x direction")
         ("H_y"     , "magnetic field in the y direction")
         ("H_z"     , "magnetic field in the z direction")
         ("H_xx"    , "nearest neighbor spin coupling Sx Sx")
         ("H_yy"    , "nearest neighbor spin exchange Sy Sy")
         ("H_zz"    , "nearest neighbor spin exchange Sz Sz")
         ("H_xx_p"  , "nearest neighbor spin coupling Sx Sx, periodic in the y direction")
         ("H_yy_p"  , "nearest neighbor spin exchange Sy Sy, periodic in the y direction")
         ("H_zz_p"  , "nearest neighbor spin exchange Sz Sz, periodic in the y direction")
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

      LatticeSite Site = SpinSite(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz"), Sp(Cell, "Sp"), Sm(Cell, "Sm");
      UnitCellOperator I(Cell, "I"); // identity operator
      InfiniteLattice Lattice(&Cell);

      UnitCellMPO XXx, YYx, ZZx, XXy, YYy, ZZy, XXyp, YYyp, ZZyp, hX, hY, hZ;
      for (int i = 0; i < Legs-1; ++i)
      {
         XXy += Sx(0)[i] * Sx(0)[i+1];
         YYy += Sy(0)[i] * Sy(0)[i+1];
         ZZy += Sz(0)[i] * Sz(0)[i+1];
      }

      for (int i = 0; i < Legs; ++i)
      {
         XXx += Sx(0)[i] * Sx(1)[i];
         YYx += Sy(0)[i] * Sy(1)[i];
         ZZx += Sz(0)[i] * Sz(1)[i];

         XXyp += Sx(0)[i] * Sx(0)[(i+1)%Legs];
         YYyp += Sy(0)[i] * Sy(0)[(i+1)%Legs];
         ZZyp += Sz(0)[i] * Sz(0)[(i+1)%Legs];

         hX += Sx(0)[i];
         hY += Sy(0)[i];
         hZ += Sz(0)[i];
      }

      Lattice["H_xx_x"] = sum_unit(XXx);
      Lattice["H_yy_x"] = sum_unit(YYx);
      Lattice["H_zz_x"] = sum_unit(ZZx);
      Lattice["H_xx_y"] = sum_unit(XXy);
      Lattice["H_yy_y"] = sum_unit(YYy);
      Lattice["H_zz_y"] = sum_unit(ZZy);
      Lattice["H_xx_yp"] = sum_unit(XXyp);
      Lattice["H_yy_yp"] = sum_unit(YYyp);
      Lattice["H_zz_yp"] = sum_unit(ZZyp);
      Lattice["H_xx"] = Lattice["H_xx_x"] + Lattice["H_xx_y"];
      Lattice["H_yy"] = Lattice["H_yy_x"] + Lattice["H_yy_y"];
      Lattice["H_zz"] = Lattice["H_zz_x"] + Lattice["H_zz_y"];
      Lattice["H_xx_p"] = Lattice["H_xx_x"] + Lattice["H_xx_yp"];
      Lattice["H_yy_p"] = Lattice["H_yy_x"] + Lattice["H_yy_yp"];
      Lattice["H_zz_p"] = Lattice["H_zz_x"] + Lattice["H_zz_yp"];
      Lattice["H_x"] = sum_unit(hX);
      Lattice["H_y"] = sum_unit(hY);
      Lattice["H_z"] = sum_unit(hZ);

      Lattice["H_J1_x"] = Lattice["H_xx_x"] + Lattice["H_yy_x"] + Lattice["H_zz_x"];
      Lattice["H_J1_y"] = Lattice["H_xx_y"] + Lattice["H_yy_y"] + Lattice["H_zz_y"];
      Lattice["H_J1_yp"] = Lattice["H_xx_yp"] + Lattice["H_yy_yp"] + Lattice["H_zz_yp"];
      Lattice["H_J1"] = Lattice["H_J1_x"] + Lattice["H_J1_y"];
      Lattice["H_J1_p"] = Lattice["H_J1_x"] + Lattice["H_J1_yp"];

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
