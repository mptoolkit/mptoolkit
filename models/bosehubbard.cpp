// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/bosehubbard.cpp
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
#include "models/boson.h"
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

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&MaxN),
          FormatDefault("Maximum number of bosons per site", MaxN).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("Bosonic 2-leg ladder with flux");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J"  , "nearest neighbor hopping")
         ("H_U"  , "on-site Coulomb repulsion N*(N-1)/2\n")
         ("H_mu" , "chemical potential = -sum_unit(N(0))")
         ("N"    , "total particle number")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = Boson(MaxN);
      UnitCell Cell(Site);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2");

      Lattice["H_J"] = sum_unit(BH(0)*B(1) + B(0)*BH(1));
      Lattice["H_U"] = sum_unit(0.5*N2(0));
      Lattice["H_mu"] = sum_unit(-N(0));
      Lattice["N"] = sum_unit(N(0));

      // information about the lattice
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
