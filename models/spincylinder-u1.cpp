// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spincylinder-u1.cpp
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
#include "models/spin-u1.h"
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
      OpDescriptions.description("U(1) spin cylinder");
      OpDescriptions.add_operators()
         ("H_J1xz"  , "nearest neighbor spin exchange z-component in the x direction")
         ("H_J1xt"  , "nearest neighbor spin exchange transverse component in the x direction")
         ("H_J1x"  , "nearest neighbor spin exchange in the x direction (H_J1xz + H_J1xt)")
         ("H_J1yz"  , "nearest neighbor spin exchange z-component in the y direction")
         ("H_J1yt"  , "nearest neighbor spin exchange transverse component in the y direction")
         ("H_J1y"  , "nearest neighbor spin exchange in the y direction (H_J1yz + H_J1yt)")
         ("H_J1"   , "nearest neighbor spin exchange H_J1x + H_J2x")
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

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(Cell);
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");

      UnitCellMPO J1xz, J1xt, J1yz, J1yt;
      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            J1xz += Sz(0)[i]*Sz(1)[i];
            J1xt += 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
            J1yz += Sz(0)[i]*Sz(0)[(i+1)%y];
            J1yt += 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%y] + Sm(0)[i]*Sp(0)[(i+1)%y]);
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            J1xz += Sz(0)[i]*Sz(0)[i+1];
            J1xt += 0.5 * (Sp(0)[i]*Sm(0)[i+1] + Sm(0)[i]*Sp(0)[i+1]);
         }
         J1xz += Sz(0)[x-1]*Sz(y+1)[0];
         J1xt += 0.5 * (Sp(0)[x-1]*Sm(y+1)[0] + Sm(0)[x-1]*Sp(y+1)[0]);
         for (int i = 0; i < x; ++i)
         {
            J1yz += Sz(0)[i]*Sz(1)[i];
            J1yt += 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
         }
      }

      Lattice["H_J1xz"] = sum_unit(J1xz);
      Lattice["H_J1xt"] = sum_unit(J1xt);
      Lattice["H_J1yz"] = sum_unit(J1yz);
      Lattice["H_J1yt"] = sum_unit(J1yt);
      Lattice["H_J1x"] = sum_unit(J1xz+J1xt);
      Lattice["H_J1y"] = sum_unit(J1yz+J1yt);
      Lattice["H_J1"] = sum_unit(J1xz+J1xt+J1yz+J1yt);

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
