// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/hubbard-u1u1.cpp
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
#include "models/fermion-u1u1.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1)xSU(2) Fermi Hubbard model");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_tup"    , "nearest neighbor hopping for up spins")
         ("H_tdown"  , "nearest neighbor hopping for down spins")
         ("H_t"      , "nearest neighbor hopping = H_tup + H_tdown")
         ("H_U"      , "on-site Coulomb interaction n_up*n_down")
         ("H_Us"     , "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
         ("H_V"      , "nearest-neighbor Coulomb interaction")
         ("H_Jup"    , "nearest-neighbor up-spin complex hopping i*(C^\\dagger_i C_{i+1} - H.c.)")
         ("H_Jdown"  , "nearest-neighbor down-spin complex hopping i*(C^\\dagger_i C_{i+1} - H.c.)")
         ("H_J"      , "nearest-neighbor complex hopping = H_Jup + H_Jdown")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      LatticeSite Site = FermionU1U1();
      UnitCell Cell(Site);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CHup(Cell, "CHup"), CHdown(Cell, "CHdown"), Cup(Cell, "Cup"),
         Cdown(Cell, "Cdown"), Pdouble(Cell, "Pdouble"),
         Hu(Cell, "Hu"), N(Cell, "N");

      Lattice["H_tup"]   = -sum_unit(dot(CHup(0), Cup(1)) - dot(Cup(0), CHup(1)));
      Lattice["H_tdown"] = -sum_unit(dot(CHdown(0), Cdown(1)) - dot(Cdown(0), CHdown(1)));
      Lattice["H_t"]     = Lattice["H_tup"] + Lattice["H_tdown"];
      Lattice["H_U"]     = sum_unit(Pdouble(0));
      Lattice["H_Us"]    = sum_unit(Hu(0));
      Lattice["H_V"]     = sum_unit(dot(N(0), N(1)));
      Lattice["H_Jup"]   = sum_unit(std::complex<double>(0,1)
                                    *(dot(CHup(0), Cup(1)) + dot(Cup(0), CHup(1))));
      Lattice["H_Jdown"] = sum_unit(std::complex<double>(0,1)
                                    *(dot(CHdown(0), Cdown(1)) + dot(Cdown(0), CHdown(1))));
      Lattice["H_J"] = Lattice["H_Jup"] + Lattice["H_Jdown"];

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
