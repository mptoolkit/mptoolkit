// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/hubbardcylinder-u1u1.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
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
#include "models/fermion-u1u1.h"
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
      OpDescriptions.set_description("U(1)xU(1) Fermi Hubbard 2D cylinder square lattice");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_txup"   , "nearest neighbor hopping in x-direction for up spins")
         ("H_txdown" , "nearest neighbor hopping in x-direction for down spins")
         ("H_tx"     , "nearest neighbor hopping in x-direction")
         ("H_tyup"   , "nearest neighbor hopping in y-direction for up spins")
         ("H_tydown" , "nearest neighbor hopping in y-direction for down spins")
         ("H_ty"     , "nearest neighbor hopping in y-direction")
         ("H_tup"    , "nearest neighbor hopping for up spins")
         ("H_tdown"  , "nearest neighbor hopping for down spins")
         ("H_t"      , "nearest neighbor hopping")
         ("H_U"      , "on-site Coulomb interaction n_up*n_down")
         ("H_Us"     , "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
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

      LatticeSite Site = FermionU1U1();
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CHup(Cell, "CHup"), CHdown(Cell, "CHdown"),
         Cup(Cell, "Cup"), Cdown(Cell, "Cdown"), Pdouble(Cell, "Pdouble"),
         Hu(Cell, "Hu");

      UnitCellMPO txup, txdown, tyup, tydown, U, Us;
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            txup += dot(CHup(0)[i], Cup(1)[i]) - dot(Cup(0)[i], CHup(1)[i]);
            txdown += dot(CHdown(0)[i], Cdown(1)[i]) - dot(Cdown(0)[i], CHdown(1)[i]);
            tyup += dot(CHup(0)[i], Cup(0)[(i+1)%y]) - dot(Cup(0)[i], CHup(0)[(i+1)%y]);
            tydown += dot(CHdown(0)[i], Cdown(0)[(i+1)%y]) - dot(Cdown(0)[i], CHdown(0)[(i+1)%y]);
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            txup += dot(CHup(0)[i], Cup(0)[i+1]) - dot(Cup(0)[i], CHup(0)[i+1]);
            txdown += dot(CHdown(0)[i], Cdown(0)[i+1]) - dot(Cdown(0)[i], CHdown(0)[i+1]);
         }
         txup += dot(CHup(0)[x-1], Cup(y+1)[0]) - dot(Cup(0)[x-1], CHup(y+1)[0]);
         txdown += dot(CHdown(0)[x-1], Cdown(y+1)[0]) - dot(Cdown(0)[x-1], CHdown(y+1)[0]);
         for (int i = 0; i < x; ++i)
         {
            tyup += dot(CHup(0)[i], Cup(1)[i]) - dot(Cup(0)[i], CHup(1)[i]);
            tydown += dot(CHdown(0)[i], Cdown(1)[i]) - dot(Cdown(0)[i], CHdown(1)[i]);
         }
      }

      for (int i = 0; i < CellSize; ++i)
      {
         U += Pdouble(0)[i];
         Us += Hu(0)[i];
      }

      Lattice["H_txup"] = -sum_unit(txup);
      Lattice["H_txdown"] = -sum_unit(txdown);
      Lattice["H_tx"] = Lattice["H_txup"] + Lattice["H_txdown"];
      Lattice["H_tyup"] = -sum_unit(tyup);
      Lattice["H_tydown"] = -sum_unit(tydown);
      Lattice["H_ty"] = Lattice["H_tyup"] + Lattice["H_tydown"];
      Lattice["H_tup"] = Lattice["H_txup"] + Lattice["H_tyup"];
      Lattice["H_tdown"] = Lattice["H_txdown"] + Lattice["H_tydown"];
      Lattice["H_t"] = Lattice["H_tup"] + Lattice["H_tdown"];
      Lattice["H_U"] = sum_unit(U);
      Lattice["H_Us"] = sum_unit(Us);

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
