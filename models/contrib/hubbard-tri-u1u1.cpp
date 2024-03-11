// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/hubbard-tri-u1u1.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
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
      OpDescriptions.set_description("U(1)xU(1) Triangular Hubbard model");
      OpDescriptions.add_operators()
         ("H_tup"    , "nearest-neighbor hopping between the apex sites of clusters for up spins")
         ("H_tdown"  , "nearest-neighbor hopping between the apex sites of clusters for down spins")
         ("H_t"      , "nearest-neighbor hopping between the apex sites of clusters")
         ("H_t2up"   , "next-nearest-neighbor hopping between the apex sites of clusters for up spins")
         ("H_t2down" , "next-nearest-neighbor hopping between the apex sites of clusters for down spins")
         ("H_t2"     , "next-nearest-neighbor hopping between the apex sites of clusters")
         ("H_tcup"   , "hopping inside the cluster for up spins")
         ("H_tcdown" , "hopping inside the cluster for down spins")
         ("H_tc"     , "hopping inside the cluster")
         ("H_tpup"   , "nearest-neighbor hopping between non-apex sites of clusters for up spins")
         ("H_tpdown" , "nearest-neighbor hopping between non-apex sites of clusters for down spins")
         ("H_tp"     , "nearest-neighbor hopping between non-apex sites of clusters")
         ("H_U"      , "on-site Coulomb interaction n_up*n_down")
         ("H_Us"     , "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
         ;
      OpDescriptions.add_cell_operators()
         ("Parity"   , "permionic swap operator for legs 0,2")
         ("R"        , "spatial reflection")
         ("N"        , "number operator")
         ("Sp"       , "spin raising operator")
         ("Sm"       , "spin lowering operator")
         ("Sz"       , "z-component of spin")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = FermionU1U1();
      UnitCell Cell(repeat(Site, 3));
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator CHup(Cell, "CHup"), CHdown(Cell, "CHdown"), Cup(Cell, "Cup"),
         Cdown(Cell, "Cdown"), Pdouble(Cell, "Pdouble"),
         Hu(Cell, "Hu"), N(Cell, "N"), Sz(Cell, "Sz"), Sp(Cell, "Sp"), Sm(Cell, "Sm"),
         R(Cell, "R");

      // parity operators
      // **NOTE** Currently the swap_gate doesn't do fermion signs, so we
      // need to add this by hand
      Cell["Parity"] = Cell.swap_gate_no_sign(0,2) * exp(math_const::pi*std::complex<double>(0,1)
                                                         *((N[0]+N[2])*N[1]+N[0]*N[2]));
      // Reflection operator.
      R = R[0]*R[1]*R[2];

      // some operators per unit cell
      N = N[0] + N[1] + N[2];
      Sz = Sz[0] + Sz[1] + Sz[2];
      Sp = Sp[0] + Sp[1] + Sp[2];
      Sm = Sm[0] + Sm[1] + Sm[2];

      Lattice["H_tup"]   = sum_unit(-(dot(CHup(0)[1], Cup(1)[1]) - dot(Cup(0)[1], CHup(1)[1])));
      Lattice["H_tdown"] = sum_unit(-(dot(CHdown(0)[1], Cdown(1)[1]) - dot(Cdown(0)[1], CHdown(1)[1])));
      Lattice["H_t"] = Lattice["H_tup"] + Lattice["H_tdown"];
      Lattice["H_t2"] = sum_unit(-(dot(CHup(0)[1], Cup(2)[1]) - dot(Cup(0)[1], CHup(2)[1])));
      Lattice["H_t2"] = sum_unit(-(dot(CHdown(0)[1], Cdown(2)[1]) - dot(Cdown(0)[1], CHdown(2)[1])));
      Lattice["H_t2"] =  Lattice["H_t2up"] + Lattice["H_t2down"];
      Lattice["H_tc"] = sum_unit(-(dot(CHup(0)[0], Cup(0)[1]) - dot(Cup(0)[0], CHup(0)[1])
                                   + dot(CHup(0)[1], Cup(0)[2]) - dot(Cup(0)[1], CHup(0)[2])
                                   + dot(CHup(0)[0], Cup(0)[2]) - dot(Cup(0)[0], CHup(0)[2])));
      Lattice["H_tc"] = sum_unit(-(dot(CHdown(0)[0], Cdown(0)[1]) - dot(Cdown(0)[0], CHdown(0)[1])
                                   + dot(CHdown(0)[1], Cdown(0)[2]) - dot(Cdown(0)[1], CHdown(0)[2])
                                   + dot(CHdown(0)[0], Cdown(0)[2]) - dot(Cdown(0)[0], CHdown(0)[2])));
      Lattice["H_tc"] =  Lattice["H_tcup"] + Lattice["H_tcdown"];
      Lattice["H_tpup"] = sum_unit(-(dot(CHup(0)[0], Cup(1)[0]) - dot(Cup(0)[0], CHup(1)[0])
                                   + dot(CHup(0)[2], Cup(1)[2]) - dot(Cup(0)[2], CHup(1)[2])));
      Lattice["H_tpdown"] = sum_unit(-(dot(CHdown(0)[0], Cdown(1)[0]) - dot(Cdown(0)[0], CHdown(1)[0])
                              + dot(CHdown(0)[2], Cdown(1)[2]) - dot(Cdown(0)[2], CHdown(1)[2])));
      Lattice["H_tp"] =  Lattice["H_tpup"] + Lattice["H_tpdown"];
      Lattice["H_U"]  = sum_unit(Pdouble(0)[0] + Pdouble(0)[1] + Pdouble(0)[2]);
      Lattice["H_Us"] = sum_unit(Hu(0)[0] + Hu(0)[1] + Hu(0)[2]);

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
