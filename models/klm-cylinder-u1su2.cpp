// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/klm-u1u1.cpp
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
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      // Parameters of the lattice (with defaults, if applicable)
      int x = 0;
      int y = 4;
      std::string FileName;

      // Options
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
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

      // Descriptions of each operator
      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1)xSU(2) Kondo Lattice Model");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_tx"   , "nearest neighbor hopping in the x direction")
         ("H_ty"   , "nearest neighbor hopping in the y direction")
         ("H_t"    , "H_tx + H_ty")
         ("H_Jx"   , "nearest-neighbor spin-spin interaction in the x direction")
         ("H_Jy"   , "nearest-neighbor spin-spin interaction in the y direction")
         ("H_J"    , "H_Jx + H_Jy")
         ("H_U"    , "on-site Coulomb interaction")
         ("H_K"    , "Kondo coupling")
         ("T_Y"    , "Translation by 1 site in the Y direction")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      // The unit cell size depends on the wrapping style. If the wrapping vector is vertical (i.e. x == 0)
      // then the unit cell size is the width y.  Otherwise the unit cell size is the wrapping offset x.
      int CellSize = x == 0 ? y : x;

      // The lattice contains only a single type of site
      LatticeSite Site = KondoU1SU2();

      // The UnitCell
      UnitCell Cell(repeat(Site, CellSize));

      // Make an infinite lattice of our unit cell
      InfiniteLattice Lattice(&Cell);

      // A short-cut to refer to an operator defined within our unit cell
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), ScSf(Cell, "ScSf"), Sf(Cell, "Sf"), Pdouble(Cell, "Pdouble"), I(Cell, "I");

      // Since the MPO's are a bit complicated, we assemble the terms on a single unit cell and then sum them
      UnitCellMPO tx, ty, Jx, Jy, U, K;

      // The xy configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            tx += dot(CH(0)[i], C(1)[i]) + dot(C(0)[i], CH(1)[i]);
            ty += dot(CH(0)[i], C(0)[(i+1)%y]) + dot(C(0)[i], CH(0)[(i+1)%y]);
            Jx += inner(Sf(0)[i], Sf(1)[i]);
            Jy += inner(Sf(0)[i], Sf(0)[(i+1)%y]);
            U += Pdouble(0)[i];
            K += ScSf(0)[i];
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            tx += dot(CH(0)[i], C(0)[i+1]) + dot(C(0)[i], CH(0)[i+1]);
            Jx += inner(Sf(0)[i], Sf(0)[i+1]);
         }
         tx += dot(CH(0)[x-1], C(y+1)[0]) + dot(C(0)[x-1], CH(y+1)[0]);
         Jx += inner(Sf(0)[x-1], Sf(y+1)[0]);
         for (int i = 0; i < x; ++i)
         {
            ty += dot(CH(0)[i], C(1)[i]) + dot(C(0)[i], CH(1)[i]);
            Jy += inner(Sf(0)[i], Sf(1)[i]);
            U += Pdouble(0)[i];
            K += ScSf(0)[i];
         }
      }

      // Define operators that have support over the infinite lattice
      Lattice["H_tx"] = -sum_unit(tx);
      Lattice["H_ty"] = -sum_unit(ty);
      Lattice["H_t"]  = Lattice["H_tx"] + Lattice["H_ty"];
      Lattice["H_Jx"] =  sum_unit(Jx);
      Lattice["H_Jy"] =  sum_unit(Jy);
      Lattice["H_J"]  = Lattice["H_Jx"] + Lattice["H_Jy"];
      Lattice["H_U"]  =  sum_unit(U);
      Lattice["H_K"]  =  sum_unit(K);

      // Translation operators.
      if (x == 0)
      {
         UnitCellMPO Trans = I(0);
         for (int i = 0; i < y-1; ++i)
         {
            Trans = Trans * Cell.swap_gate(i, i+1);
         }
         Lattice["Ty"] = prod_unit_left_to_right(Trans, CellSize);
      }

      if (x == 0)
      {
         UnitCellMPO Trans = I(0);
         for (int i = 0; i < y-1; ++i)
         {
            Trans = Trans * Cell.swap_gate_no_sign(i, i+1);
         }
         Lattice["sTy"] = prod_unit_left_to_right(Trans, CellSize);
      }


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
