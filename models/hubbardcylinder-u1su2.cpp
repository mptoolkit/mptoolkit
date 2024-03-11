// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/hubbardcylinder-u1su2.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
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
#include "models/fermion-u1su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      int x = 0;
      int y = 4;

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

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1)xSU(2) Fermi Hubbard 2D cylinder square lattice");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_tx" , "nearest neighbor hopping in y-direction")
         ("H_ty" , "nearest neighbor hopping in x-direction")
         ("H_t"  , "nearest neighbor hopping")
         ("H_U"  , "on-site Coulomb interaction n_up*n_down")
         ("H_Us" , "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      int CellSize = x == 0 ? y : x;

      LatticeSite Site = FermionU1SU2();
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), Pdouble(Cell, "Pdouble"),
         Hu(Cell, "Hu"), N(Cell, "N");

      UnitCellMPO tx, ty, Pd, hu;
      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            tx += dot(CH(0)[i], C(1)[i]) + dot(C(0)[i], CH(1)[i]);
            ty += dot(CH(0)[i], C(0)[(i+1)%y]) + dot(C(0)[i], CH(0)[(i+1)%y]);
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            tx += dot(CH(0)[i], C(0)[i+1]) + dot(C(0)[i], CH(0)[i+1]);
         }
         tx += dot(CH(0)[x-1], C(y+1)[0]) + dot(C(0)[x-1], CH(y+1)[0]);
         for (int i = 0; i < x; ++i)
         {
            ty += dot(CH(0)[i], C(1)[i]) + dot(C(0)[i], CH(1)[i]);
         }
      }

      for (int i = 0; i < CellSize; ++i)
      {
         Pd += Pdouble(0)[i];
         hu += Hu(0)[i];
      }

      Lattice["H_tx"] = sum_unit(tx);
      Lattice["H_ty"] = sum_unit(ty);
      Lattice["H_t"] = sum_unit(tx+ty);
      Lattice["H_Us"] = sum_unit(Pd);
      Lattice["H_U"] = sum_unit(hu);

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
