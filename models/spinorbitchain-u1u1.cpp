// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinchain-u1.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1)xU(1) S-T spin chain");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_Sz"  , "nearest neighbor spin coupling Sz Sz")
         ("H_St"  , "nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
         ("H_S"   , "nearest neighbor spin = H_Sz + H_St")
         ("H_Tz"  , "nearest neighbor pseudospin coupling Tz Tz")
         ("H_Tt"  , "nearest neighbor pseudospin exchange (1/2)(Tp Tm + Tm Tp)")
         ("H_T"   , "nearest neighbor pseudospin = H_Tz + H_Tt")
         ("H_ST"  , "nearest neighbor spin * pseudospin")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      LatticeSite SiteS = SpinU1(Spin, "Sz");
      LatticeSite SiteT = SpinU1(Spin, "Tz");
      SymmetryList ST("Sz:U(1),Tz:U(1)");
      UnitCell Cell(ST, SiteS, SiteT);

      // On the 2-site unit cell, S^a[0] is the S sites, S^a[1] is the T sites
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      UnitCellOperator Tp(Cell, "Tp"), Tm(Cell, "Tm"), Tz(Cell, "Tz");

      // Unit cell operators S and T
      Sp(0) = Sp(0)[0];
      Sm(0) = Sm(0)[0];
      Sz(0) = Sz(0)[0];

      Tp(0) = Sp(0)[1];
      Tm(0) = Sm(0)[1];
      Tz(0) = Sz(0)[1];

      InfiniteLattice Lattice(&Cell);

      Lattice["H_Sz"] = sum_unit(Sz(0)*Sz(1));
      Lattice["H_St"] = 0.5 * sum_unit(Sp(0)*Sm(1) + Sm(0)*Sp(1));
      Lattice["H_S"]  = Lattice["H_Sz"] + Lattice["H_St"];

      Lattice["H_Tz"] = sum_unit(Tz(0)*Tz(1));
      Lattice["H_Tt"] = 0.5 * sum_unit(Tp(0)*Tm(1) + Tm(0)*Tp(1));
      Lattice["H_T"]  = Lattice["H_Tz"] + Lattice["H_Tt"];

      Lattice["H_ST"] = sum_unit( ( Sz(0)*Sz(1) + 0.5 * (Sp(0)*Sm(1) + Sm(0)*Sp(1)) )
                                 *( Tz(0)*Tz(1) + 0.5 * (Tp(0)*Tm(1) + Tm(0)*Tp(1)) ) );

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
