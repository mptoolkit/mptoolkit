// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/z2-lgt-u1.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2020-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1) 1+1D Z2 lattice gauge theory

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"
#include "models/spinlessantifermion-u1.h"
#include "models/spin.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      half_int Spin = 0.5;
      bool Tau = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         //("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the link spins", Spin).c_str()) // This should only ever be 0.5.
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) 1+1D Z2 lattice gauge theory");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("G"    , "Gauss's law operator")
         ;
      OpDescriptions.add_operators()
         ("H_t"  , "nearest-neighbor hopping")
         ("H_m"  , "fermion mass")
         ("H_x"  , "magnetic field in the x direction")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         std::cerr << "Operators:" << std::endl << OpDescriptions;
         return 1;
      }

      LatticeSite FSite = SpinlessFermionU1();
      LatticeSite SSite = SpinSite(Spin);

      UnitCell Cell(FSite.GetSymmetryList(), FSite, SSite);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N"), P(Cell, "P"), I(Cell, "I"),
                       X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");

      UnitCellMPO t, m, x_field;

      t += Z(0)[1] * dot(CH(0)[0], C(1)[0]) + Z(0)[1] * dot(CH(1)[0], C(0)[0]);
      m += N(0)[0];
      x_field += X(0)[1];

      Lattice["H_t"] = sum_unit(t);
      Lattice["H_m"] = sum_unit(m);
      Lattice["H_x"] = sum_unit(x_field);

      // Gauss's law operators.
      UnitCellOperator G(Cell, "G");

      G = P(0)[0] * X(0)[1] * X(-1)[1];

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
