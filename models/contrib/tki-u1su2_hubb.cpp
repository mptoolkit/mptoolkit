// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/tki-u1su2_hubb.cpp
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
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
#include "models/spin-su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"
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
    OpDescriptions.add_operators()
      ("H_t"  , "nearest neighbour fermion hopping")
      ("H_U"  , "on-site Coulomb interaction n_up*n_down")
      ("H_J1" , "nearest neighbour spin exchange")
      ("H_K"  , "Kondo coupling between fermion and spin")
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
    LatticeSite fSite = SpinSU2(0.5);
    UnitCell Cell(cSite.GetSymmetryList(), cSite, fSite);
    InfiniteLattice Lattice(&Cell);

    UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), S(Cell, "S"),
       p(Cell, "p"), pH(Cell, "pH"), Pi(Cell, "Pi"), Pdouble(Cell, "Pdouble");

    // note: need to define this *BEFORE* constructing the InfiniteLattice object
    p(0) = std::sqrt(0.5) * (C(1)[0] - C(-1)[0]);
    pH(0) = adjoint(p(0));

    Pi(0) = outer(pH(0), p(0));

    Lattice["H_t"]  = sum_unit(dot(CH(0)[0], C(1)[0]) + dot(C(0)[0], CH(1)[0]));
    Lattice["H_U"]  = sum_unit(Pdouble(0)[0]);
    Lattice["H_J1"] = sum_unit(inner(S(0)[1], S(1)[1]));
    Lattice["H_K"]  = sum_unit(inner(S(0)[1], Pi(0)));

    // Information about the lattice
    Lattice.set_description("U(1)xSU(2) Kondo lattice model");
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
