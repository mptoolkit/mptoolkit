// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/bosehubbard-u1.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "models/boson-u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      int MaxN = DefaultMaxN;
      bool QLM = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&MaxN),
          FormatDefault("Maximum number of bosons per site", MaxN).c_str())
         ("qlm", prog_opt::bool_switch(&QLM),
          "include terms for the mapping of the QLM")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) Bose-Hubbard model");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J"    , "nearest-neighbor hopping")
         ("H_Jm"   , "nearest-neighbor hopping in negative direction")
         ("H_Jp"   , "nearest-neighbor hopping in positive direction")
         ("H_U"    , "on-site Coulomb repulsion N*(N-1)/2")
         ("H_V"    , "nearest-neighbour Coulomb repulsion")
         ("H_delta", "QLM gauge site potential", "QLM enabled",
         [&QLM]()->bool{return QLM;})
         ("H_W"    , "next-nearest-neighbour gauge site Coulomb repulsion", "QLM enabled",
         [&QLM]()->bool{return QLM;})
         ("H_chi"  , "QLM staggering potential", "QLM enabled",
         [&QLM]()->bool{return QLM;})
         ;
      OpDescriptions.add_functions()
         ("H_J"    , "nearest-neighbor complex hopping")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = BosonU1(MaxN);
      UnitCell Cell(Site);
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2");

      InfiniteLattice Lattice(&Cell);

      Lattice["H_J"] = sum_unit(BH(0)*B(1) + B(0)*BH(1));

      Lattice["H_Jm"] = sum_unit(BH(0)*B(1));
      Lattice["H_Jp"] = sum_unit(B(0)*BH(1));
      Lattice.func("H_J")(arg("theta")) = "exp(-i*theta)*H_Jm + exp(i*theta)*H_Jp";

      Lattice["H_U"] = sum_unit(0.5*N2(0));
      Lattice["H_V"] = sum_unit(N(0)*N(1));

      if (QLM)
      {
         Lattice["H_delta"] = sum_unit(N(1), 2);
         Lattice["H_W"] = sum_unit(N(1)*N(3), 2);
         Lattice["H_chi"] = sum_unit(0.5*(N(1)-N(3)), 4);
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
