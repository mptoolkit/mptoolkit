// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/hubbardholstein-u1su2.cpp
//
// Copyright (C) 2015-2024 Ian McCulloch <ian@qusim.net>
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
#include "models/fermion-u1su2.h"
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
      int MaxN = 4;

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
      OpDescriptions.description("U(1)xSU(2) Hubbard-Holstein model");
      OpDescriptions.author("Ian McCulloch", "ian@qusim.net");
      OpDescriptions.add_operators()
         ("H_t"   , "nearest neighbour fermion hopping")
         ("H_U"   , "on-site coulomb repulsion for the fermions")
         ("H_Us"  , "on-site coulomb repulsion for the fermions, symmetric version")
         ("H_w"   , "on-site phonon energy")
         ("H_g"   , "Fermion-phonon coupling")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      LatticeSite cSite = FermionU1SU2();
      LatticeSite pSite = Boson(MaxN);
      UnitCell Cell(cSite.GetSymmetryList(), cSite, pSite);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), S(Cell, "S"), BH(Cell, "BH"), B(Cell, "B"),
                       Pdouble(Cell, "Pdouble"), Hu(Cell, "Hu"), N(Cell, "N"), I(Cell, "I");

      // FIXME: the hopping term should have a negative sign.
      // This is not important for a bipartite lattice, but will make a difference
      // if a NNN hopping term is included.
      Lattice["H_t"]  = sum_unit(dot(CH(0)[0], C(1)[0]) + adjoint(dot(CH(0)[0], C(1)[0])));
      Lattice["H_U"]  =  sum_unit(Pdouble(0)[0]);
      Lattice["H_Us"] =  sum_unit(Hu(0)[0]);
      Lattice["H_w"]  = sum_unit(N(0)[1]); //dot(BH(0)[1], B(0)[1]));
      Lattice["H_g"]  = sum_unit((N(0)[0] - I(0)[0]) * (BH(0)[1] + B(0)[1]));
      Lattice["H_g0"]  = sum_unit(N(0)[0] * (BH(0)[1] + B(0)[1]));

      // Information about the lattice
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice to disc
      pheap::ExportObject(FileName, Lattice);
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
