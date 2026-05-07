// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spinladder-u1.cpp
//
// Copyright (C) 2004-2026 Ian McCulloch <ian@qusim.net>
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
#include "models/spin-u1.h"
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
      OpDescriptions.description("U(1) spin ladder");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1z_x" , "nearest neighbor spin exchange z-component in the x (infinite) direction")
         ("H_J1t_x" , "nearest neighbor spin exchange transverse component in the x (infinite) direction")
         ("H_J1_x"  , "nearest neighbor spin exchange in the x (infinite) direction")
         ("H_J1z_y" , "nearest neighbor spin exchange z-component in the y (rung) direction")
         ("H_J1t_y" , "nearest neighbor spin exchange transverse component in the y (rung) direction")
         ("H_J1_y"  , "nearest neighbor spin exchange in the y (rung) direction")
         ("H_J1z_yp", "nearest neighbor spin exchange z-component in the y (rung) direction, including periodic term")
         ("H_J1t_yp", "nearest neighbor spin exchange transverse component in the y (rung) direction, including periodic term")
         ("H_J1_yp" , "nearest neighbor spin exchange in the y (rung) direction, including periodic term")
         ("H_J1z"   , "nearest neighbor spin exchange z-component H_J1z_x + H_J1z_y")
         ("H_J1t"   , "nearest neighbor spin exchange transverse component H_J1t_x + H_J1t_y")
         ("H_J1"    , "nearest neighbor spin exchange H_J1_x + H_J1_y")
         ("H_J1z_p" , "nearest neighbor spin exchange z-component H_J1z_x + H_J1z_yp")
         ("H_J1t_p" , "nearest neighbor spin exchange transverse component H_J1t_x + H_J1t_yp")
         ("H_J1_p"  , "nearest neighbor spin exchange H_J1_x + H_J1_yp (periodic in y direction)")
         ("H_J1x"   , "alias for H_J1_x")
         ("H_J1y"   , "alias for H_J1_y")
         ("H_J1yp"  , "alias for H_J1_yp")
         ("H_J1p"   , "alias for H_J1_p")
         ("H_z"     , "magnetic field in the z direction")
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

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      InfiniteLattice Lattice(&Cell);

      UnitCellMPO J1zx, J1tx, J1zy, J1ty, J1zyp, J1typ, hZ;
      for (int i = 0; i < Legs-1; ++i)
      {
         J1zy += Sz(0)[i] * Sz(0)[i+1];
         J1ty += 0.5 * (Sp(0)[i] * Sm(0)[i+1] + Sm(0)[i] * Sp(0)[i+1]);
      }

      for (int i = 0; i < Legs; ++i)
      {
         J1zx += Sz(0)[i] * Sz(1)[i];
         J1tx += 0.5 * (Sp(0)[i] * Sm(1)[i] + Sm(0)[i] * Sp(1)[i]);

         J1zyp += Sz(0)[i] * Sz(0)[(i+1)%Legs];
         J1typ += 0.5 * (Sp(0)[i] * Sm(0)[(i+1)%Legs] + Sm(0)[i] * Sp(0)[(i+1)%Legs]);

         hZ += Sz(0)[i];
      }

      Lattice["H_J1z_x"] = sum_unit(J1zx);
      Lattice["H_J1t_x"] = sum_unit(J1tx);
      Lattice["H_J1_x"] = Lattice["H_J1z_x"] + Lattice["H_J1t_x"];
      Lattice["H_J1z_y"] = sum_unit(J1zy);
      Lattice["H_J1t_y"] = sum_unit(J1ty);
      Lattice["H_J1_y"] = Lattice["H_J1z_y"] + Lattice["H_J1t_y"];
      Lattice["H_J1z_yp"] = sum_unit(J1zyp);
      Lattice["H_J1t_yp"] = sum_unit(J1typ);
      Lattice["H_J1_yp"] = Lattice["H_J1z_yp"] + Lattice["H_J1t_yp"];
      Lattice["H_J1z"] = Lattice["H_J1z_x"] + Lattice["H_J1z_y"];
      Lattice["H_J1t"] = Lattice["H_J1t_x"] + Lattice["H_J1t_y"];
      Lattice["H_J1"] = Lattice["H_J1_x"] + Lattice["H_J1_y"];
      Lattice["H_J1z_p"] = Lattice["H_J1z_x"] + Lattice["H_J1z_yp"];
      Lattice["H_J1t_p"] = Lattice["H_J1t_x"] + Lattice["H_J1t_yp"];
      Lattice["H_J1_p"] = Lattice["H_J1_x"] + Lattice["H_J1_yp"];
      Lattice["H_J1x"] = Lattice["H_J1_x"];
      Lattice["H_J1y"] = Lattice["H_J1_y"];
      Lattice["H_J1yp"] = Lattice["H_J1_yp"];
      Lattice["H_J1p"] = Lattice["H_J1_p"];
      Lattice["H_z"] = sum_unit(hZ);

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
