// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/bosehubbard-ladder-u1.cpp
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
      int MaxN = 5;
      int Width = 2;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&MaxN),
          FormatDefault("Maximum number of bosons per site", MaxN).c_str())
         ("width,w", prog_opt::value(&Width),
          FormatDefault("Width of the ladder", Width).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) Bose-Hubbard ladder");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J"    , "nearest-neighbor hopping on the legs")
         ("H_K"    , "nearest-neighbor hopping on the rungs")
         ("H_U"    , "on-site Coulomb repulsion N*(N-1)/2")
         ("H_U12"  , "nearest-neighbor Coulomb repulsion on the rungs")
         ("D"      , "difference in occupation number between leg (W-1) and leg 0")
         ("D2"     , "squared difference in occuptation number between leg (W-1) and leg 0")
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
      UnitCell Cell(repeat(Site, Width));
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2"),
         Delta(Cell, "D");

      // difference of boson occupation numbers between edges of the ladder
      Delta(0) = N(0)[Width-1] - N(0)[0];
      Delta.set_description("difference in occupation number between leg (W-1) and leg 0");


      UnitCellMPO HJ, HK, HU, HU12;
      for (int i = 0; i < Width; ++i)
      {
         HJ -= BH(0)[i]*B(1)[i] + B(0)[i]*BH(1)[i];
         HU += 0.5*N2(0)[i];
      }

      for (int i = 0; i < Width-1; ++i)
      {
         HK -= BH(0)[i]*B(0)[i+1] + B(0)[i]*BH(0)[i+1];
         HU12 = N(0)[i] * N(0)[i+1];
      }

      Lattice["H_J"]   = sum_unit(HJ);
      Lattice["H_K"]   = sum_unit(HK);
      Lattice["H_U"]   = sum_unit(HU);
      Lattice["H_U12"] = sum_unit(HU12);
      Lattice["D"]     = sum_unit(Delta);
      Lattice["D2"]    = sum_unit(Delta*Delta);

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
