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
    int w;
    std::string FileName;

    prog_opt::options_description desc("Allowed options", terminal::columns());
    desc.add_options()
      ("help", "show this help message")
      ("width,w", prog_opt::value(&w), "cylinder width [must be >= 3]")
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
      ("H_K"  , "Diagonal Kondo coupling between fermion and spin")
      ("H_K1" , "Perpendicular Kondo coupling between fermion and spin")
      ;
    OpDescriptions.add_cell_operators()
       ("p"   , "p-wave annihilation")
       ("pH"  , "p-wave creation")
       ("Pi"  , "p-wave spin vector")
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
//    UnitCell Cell(cSite.GetSymmetryList(), cSite, fSite);
    UnitCell Site (cSite.GetSymmetryList(), cSite, fSite);
    UnitCell Cell(repeat(Site, w));
    InfiniteLattice Lattice(&Cell);
    UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), S(Cell, "S"), 
       p(Cell, "p"), pH(Cell, "pH"), Pi(Cell, "Pi"), Pdouble(Cell, "Pdouble");

    // note: need to define this *BEFORE* constructing the InfiniteLattice object
    UnitCellMPO H_t_cell, H_U_cell, H_J1_cell, H_K_cell, H_K1_cell;
    // Intrachain Hamiltonians
    for (int i = 0; i <= 2*w-1; i+=2)
    {
      H_t_cell += dot(CH(0)[i], C(1)[i]) + dot(C(0)[i], CH(1)[i]);
      H_U_cell += Pdouble(0)[i];
      H_J1_cell += inner(S(0)[i+1], S(1)[i+1]);

      p(i) = std::sqrt(0.5) * (C(1)[i] - C(-1)[i]);
      pH(i) = adjoint(p(i));
      Pi(i) = outer(pH(i), p(i));
      H_K_cell += inner(S(0)[i+1], Pi(i));

      H_K1_cell += inner(S(0)[i], S(0)[i+1]);
    }

    // Interchain coupling Hamiltonians
    for (int i = 3; i < 2*w-1; i+=2)
    {
      p(i) = std::sqrt(0.5) * (C(0)[i+1] - C(0)[i-3]);
      pH(i) = adjoint(p(i));
      Pi(i) = outer(pH(i), p(i));
      H_K_cell += inner(S(0)[i], Pi(i));
    }

    for (int i = 0; i < 2*(w-1); i+=2)
    {
      H_t_cell += dot(CH(0)[i], C(0)[i+2]) + dot(C(0)[i], CH(0)[i+2]);
      H_J1_cell += inner(S(0)[i+1], S(0)[i+3]);
    }
    // Couple last Hubbard chain to first Hubbard chain
    H_t_cell += dot(CH(0)[2*(w-1)], C(0)[0]) + dot(C(0)[2*(w-1)], CH(0)[0]);
    // Couple last Heisenberg chain to first Heisenberg chain
    H_J1_cell += inner(S(0)[2*w-1], S(0)[1]);

    // Couple last set of chains to second last and first set of chains
    p(0) = std::sqrt(0.5) * (C(0)[0] - C(0)[2*w-4]);
    pH(0) = adjoint(p(0));
    Pi(0) = outer(pH(0), p(0));
    H_K_cell += inner(S(0)[2*w-1], Pi(0));

    // Couple first set of chains to second set of last set ofchains
    p(0) = std::sqrt(0.5) * (C(0)[2] - C(0)[2*(w-1)]);
    pH(0) = adjoint(p(0));
    Pi(0) = outer(pH(0), p(0));
    H_K_cell += inner(S(0)[1], Pi(0));

    Lattice["H_t"]  = sum_unit(H_t_cell);
    Lattice["H_U"]  = sum_unit(H_U_cell);
    Lattice["H_J1"] = sum_unit(H_J1_cell);
    Lattice["H_K"]  = sum_unit(H_K_cell);
    Lattice["H_K1"] = sum_unit(H_K1_cell);

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
