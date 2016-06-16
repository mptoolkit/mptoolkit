// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinchain-su2.cpp
//
// Copyright (C) 2016 Jason Pillay <pillayjason@hotmail.com>
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
#include "models/spin-u1.h"
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
    OpDescriptions.description("U(1)xU(1) p-wave Kondo lattice model");
    OpDescriptions.author("Jason Pillay", "pillayjason@hotmail.com");
    OpDescriptions.add_operators()
      ("H_t"  , "nearest neighbour fermion hopping")
      ("H_J1" , "nearest neighbour spin exchange")
      ("H_K"  , "Kondo coupling between fermion and spin")
      ;
    OpDescriptions.add_cell_operators()
       ("pup"   , "annihilation up spin p-wave")
       ("pdown" , "annihilation down spin p-wave")
       ("pHup"  , "creation up spin p-wave")
       ("pHdown", "creation down spin p-wave")
       ("Pi_z"  , "p-wave z-component")
       ("Pi_p" , "p-wave raising operator")
       ("Pi_m" , "p-wave lowering operator")
       ;

    if (vm.count("help") || !vm.count("out"))
    {
      print_copyright(std::cerr);
      std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
      std::cerr << desc << '\n';
      std::cerr << OpDescriptions;
      return 1;
    }

    LatticeSite cSite = FermionU1U1();
    LatticeSite fSite = SpinU1(0.5);
    UnitCell Cell(cSite.GetSymmetryList(), cSite, fSite);
    UnitCellOperator CHup(Cell, "CHup"), CHdown(Cell, "CHdown"), Cup(Cell, "Cup"), 
       Cdown(Cell, "Cdown"), Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz"),
       pup(Cell, "pup"), pdown(Cell, "pdown"), pHup(Cell, "pHup"), pHdown(Cell, "pHdown"),
       Pi_z(Cell, "Pi_z"), Pi_p(Cell, "Pi_p"), Pi_m(Cell, "Pi_m");

    // note: need to define this *BEFORE* constructing the InfiniteLattice object
    pup(0) = std::sqrt(0.5) * (Cup(1)[0] - Cup(-1)[0]);
    pdown(0) = std::sqrt(0.5) * (Cdown(1)[0] - Cdown(-1)[0]);
    pHup(0) = adjoint(pup(0));
    pHdown(0) = adjoint(pdown(0));

    Pi_z(0) = 0.5 * (pHup(0)*pup(0) - pHdown(0)*pdown(0));
    Pi_p(0) = pHup(0)*pdown(0);
    Pi_m(0) = pHdown(0)*pup(0);

    InfiniteLattice Lattice(Cell);

    Lattice["H_t"] = sum_unit(CHup(0)[0]*Cup(1)[0] - Cup(0)[0]*CHup(1)[0]
		   + CHdown(0)[0]*Cdown(1)[0] - Cdown(0)[0]*CHdown(1)[0]);
    Lattice["H_J1"] = sum_unit(Sz(0)[1]*Sz(1)[1] + 0.5*(Sp(0)[1]*Sm(1)[1] + Sm(0)[1]*Sp(1)[1]));
    Lattice["H_K"] = 0.5*(sum_unit(Sp(0)[1] * Pi_m(0)) + sum_unit(Sm(0)[1] * Pi_p(0)))
		   + sum_unit(Sz(0)[1] * Pi_z(0));

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


