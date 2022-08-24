// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinchain-u1u1u1-subset-su4.cpp
//
// Copyright (C) 2004-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "models/contrib/spin-u1u1u1-subset-su4.h"
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
      OpDescriptions.set_description("SU(4) spin chain in the U(1)xU(1)xU(1) basis");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_t"  , "nearest neighbor T-channel exchange")
         ("H_v"  , "nearest neighbor V-channel exchange")
         ("H_u"  , "nearest neighbor U-channel exchange")
         ("H_w"  , "nearest neighbor T-channel exchange")
         ("H_x"  , "nearest neighbor V-channel exchange")
         ("H_y"  , "nearest neighbor U-channel exchange")
         ("H_3"  , "nearest neighbor L3-L3")
         ("H_8"  , "nearest neighbor L8-L8")
         ("H_15" , "nearest neighbor L15-L15")
         ("H_SU4", "0.25*(sum over SU(4) generators) = 0.25*(H_t + H_v + H_u + H_w + H_x + H_y + H_3 + H_8 + H_15)")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      LatticeSite Site = SpinU1U1U1();
      UnitCell Cell(Site);
      UnitCellOperator Tp(Cell, "Tp"), Tm(Cell, "Tm"), Vp(Cell, "Vp"), Vm(Cell, "Vm"), Up(Cell, "Up"), Um(Cell, "Um");
      UnitCellOperator Wp(Cell, "Wp"), Wm(Cell, "Wm"), Xp(Cell, "Xp"), Xm(Cell, "Xm"), Yp(Cell, "Yp"), Ym(Cell, "Ym");
      UnitCellOperator L3(Cell, "L3"), L8(Cell, "L8"), L15(Cell, "L15"), I(Cell, "I");
      InfiniteLattice Lattice(&Cell);

      Lattice["H_t"] = 0.5 * sum_unit(Tp(0)*Tm(1) + Tm(0)*Tp(1));
      Lattice["H_v"] = 0.5 * sum_unit(Vp(0)*Vm(1) + Vm(0)*Vp(1));
      Lattice["H_u"] = 0.5 * sum_unit(Up(0)*Um(1) + Um(0)*Up(1));
      Lattice["H_w"] = 0.5 * sum_unit(Wp(0)*Wm(1) + Wm(0)*Wp(1));
      Lattice["H_x"] = 0.5 * sum_unit(Xp(0)*Xm(1) + Xm(0)*Xp(1));
      Lattice["H_y"] = 0.5 * sum_unit(Yp(0)*Ym(1) + Ym(0)*Yp(1));
      Lattice["H_3"] = sum_unit(L3(0)*L3(1));
      Lattice["H_8"] = sum_unit(L8(0)*L8(1));
      Lattice["H_15"] = sum_unit(L15(0)*L15(1));
      Lattice["H_SU4"] = 0.25*(Lattice["H_t"] + Lattice["H_v"] + Lattice["H_u"] +
                                 Lattice["H_w"] + Lattice["H_x"] + Lattice["H_y"] +
                                 Lattice["H_3"] + Lattice["H_8"] + Lattice["H_15"]);

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
