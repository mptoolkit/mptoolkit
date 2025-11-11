// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/tasaki-u1u1.cpp
//
// Copyright (C) 2025 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1)xU(1) 1D Tasaki model

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/fermion-u1u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      double nu = 1.0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("nu", prog_opt::value(&nu), FormatDefault("nu", nu).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1)xU(1) Tasaki model");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_t" , "nearest-neighbor hopping")
         ("H_U" , "on-site Coulomb interaction n_up*n_down")
         ("H_Us", "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         std::cerr << "Operators:" << std::endl << OpDescriptions;
         return 1;
      }

      LatticeSite Site = FermionU1U1();
      UnitCell Cell(Site.GetSymmetryList(), Site, Site);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CHup(Cell, "CHup"), CHdown(Cell, "CHdown"), Cup(Cell, "Cup"),
                       Cdown(Cell, "Cdown"), Pdouble(Cell, "Pdouble"),
                       Hu(Cell, "Hu"), N(Cell, "N"), I(Cell, "I");

      UnitCellMPO AHup   = CHup(0)[0]   - nu * (CHup(0)[1] + CHup(-1)[1]);
      UnitCellMPO AHdown = CHdown(0)[0] - nu * (CHdown(0)[1] + CHdown(-1)[1]);
      UnitCellMPO Aup    = Cup(0)[0]    - nu * (Cup(0)[1] + Cup(-1)[1]);
      UnitCellMPO Adown  = Cdown(0)[0]  - nu * (Cdown(0)[1] + Cdown(-1)[1]);

      AHup.set_description("AHup");
      Lattice.GetUnitCell().assign_operator("AHup", AHup);
      AHdown.set_description("AHdown");
      Lattice.GetUnitCell().assign_operator("AHdown", AHdown);
      Aup.set_description("Aup");
      Lattice.GetUnitCell().assign_operator("Aup", Aup);
      Adown.set_description("Adown");
      Lattice.GetUnitCell().assign_operator("Adown", Adown);

      UnitCellMPO BHup   = CHup(0)[1]   + nu * (CHup(0)[0] + CHup(1)[0]);
      UnitCellMPO BHdown = CHdown(0)[1] + nu * (CHdown(0)[0] + CHdown(1)[0]);
      UnitCellMPO Bup    = Cup(0)[1]    + nu * (Cup(0)[0] + Cup(1)[0]);
      UnitCellMPO Bdown  = Cdown(0)[1]  + nu * (Cdown(0)[0] + Cdown(1)[0]);

      BHup.set_description("BHup");
      Lattice.GetUnitCell().assign_operator("BHup", BHup);
      BHdown.set_description("BHdown");
      Lattice.GetUnitCell().assign_operator("BHdown", BHdown);
      Bup.set_description("Bup");
      Lattice.GetUnitCell().assign_operator("Bup", Bup);
      Bdown.set_description("Bdown");
      Lattice.GetUnitCell().assign_operator("Bdown", Bdown);

      UnitCellMPO t = dot(BHup, Bup) + dot(BHdown, Bdown);
      UnitCellMPO U = Pdouble(0)[0] + Pdouble(0)[1];
      UnitCellMPO Us = Hu(0)[0] + Hu(0)[1];

      // FIXME: the hopping term should have a negative sign.
      // This is not important for a bipartite lattice, but will make a difference
      // if a NNN hopping term is included.
      Lattice["H_t"] = sum_unit(t);
      Lattice["H_U"] = sum_unit(U);
      Lattice["H_Us"] = sum_unit(Us);

      // Information about the lattice
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice to disk
      pheap::ExportObject(FileName, Lattice);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      return 1;
   }
}
