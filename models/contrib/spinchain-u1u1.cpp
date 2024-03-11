// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spinchain-u1u1.cpp
//
// Copyright (C) 2004-2022 Ian McCulloch <ian@qusim.net>
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
#include "models/contrib/spin-u1u1.h"
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
      OpDescriptions.set_description("SU(3) Lai-Sutherland Spin-1 chain in the U(1)xU(1) basis");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_t"  , "nearest neighbor T-channel exchange")
         ("H_v"  , "nearest neighbor V-channel exchange")
         ("H_u"  , "nearest neighbor U-channel exchange")
         ("H_3"  , "nearest neighbor L3-L3")
         ("H_8"  , "nearest neighbor L8-L8")
         ("H_BQ" , "Bilinear-biquadratic model = (4/3) + 0.5*(H_t + H_v + H_u + H_3 + H_8)")
         ("H_SU3" , "Bilinear-biquadratic model = 0.25*(H_t + H_v + H_u + H_3 + H_8) = 0.25*(sum over SU(3) generators)")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      LatticeSite Site = SpinU1U1();
      UnitCell Cell(Site);
      UnitCellOperator Tp(Cell, "Tp"), Tm(Cell, "Tm"), Vp(Cell, "Vp"), Vm(Cell, "Vm"), Up(Cell, "Up"), Um(Cell, "Um");
      UnitCellOperator Sz(Cell, "Sz"), Qz(Cell, "Qz"), L3(Cell, "L3"), L8(Cell, "L8"), I(Cell, "I");
      InfiniteLattice Lattice(&Cell);

      Lattice["H_t"] = 0.5 * sum_unit(Tp(0)*Tm(1) + Tm(0)*Tp(1));
      Lattice["H_v"] = 0.5 * sum_unit(Vp(0)*Vm(1) + Vm(0)*Vp(1));
      Lattice["H_u"] = 0.5 * sum_unit(Up(0)*Um(1) + Um(0)*Up(1));
      Lattice["H_3"] = sum_unit(L3(0)*L3(1));
      Lattice["H_8"] = sum_unit(L8(0)*L8(1));
      Lattice["H_BQ"] = sum_unit(I(0)*(4.0/3.0)) + 0.5*(Lattice["H_t"] + Lattice["H_v"] + Lattice["H_u"] + Lattice["H_3"] + Lattice["H_8"]);
      Lattice["H_SU3"] = 0.25*(Lattice["H_t"] + Lattice["H_v"] + Lattice["H_u"] + Lattice["H_3"] + Lattice["H_8"]);

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
