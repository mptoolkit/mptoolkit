// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/z2-lgt-cylinder-u1.cpp
//
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1) 2+1D Z2 lattice guage theory
// Example for (x,y) = (0,2).
//
//     UC 1          UC 2
// 0--1--6--7--0--1--
// |     |     |
// 5     11    5
// |     |     |
// 3--4--9--10-3--4--  ...
// |     |     |
// 2     8     2
// |     |     |
// 0--1--6--7--0--1--
//
// 0 and 9 are fermions, while 3 and 6 are antifermions.
// The other sites are link sites represented by spin-1/2s.

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
      int x = 0;
      int y = 4;
      half_int Spin = 0.5;
      bool Bosonic = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         //(",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str()) TODO: twisted boundaries
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector (must be even)", y).c_str())
         //("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the link spins", Spin).c_str()) // This should only ever be 0.5.
         ("bosonic", prog_opt::bool_switch(&Bosonic), "use hardcore bosons instead of spinless fermions on the matter sites")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) 2+1D Z2 lattice gauge theory");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_tx" , "nearest-neighbor hopping in y-direction")
         ("H_ty" , "nearest-neighbor hopping in x-direction")
         ("H_t"  , "nearest-neighbor hopping")
         ("H_m"  , "fermion mass")
         ("H_J"  , "plaquette interactions")
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

      int CellSize = x == 0 ? 3*y : 3*x;

      LatticeSite FSite = SpinlessFermionU1("N", "P", Bosonic);
      LatticeSite SSite = SpinSite(Spin);
      UnitCell UCell(FSite.GetSymmetryList(), FSite, SSite, SSite);
      UnitCell Cell(repeat(UCell, CellSize/3));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N"), P(Cell, "P"), I(Cell, "I"),
                       X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");

      UnitCellMPO tx, ty, m, J, x_field;
      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < 3*y; i += 3)
         {
            tx += Z(0)[i+1] * dot(CH(0)[i], C(1)[i])           + Z(0)[i+1] * dot(CH(1)[i],           C(0)[i]);
            ty += Z(0)[i+2] * dot(CH(0)[i], C(0)[(i+3)%(3*y)]) + Z(0)[i+2] * dot(CH(0)[(i+3)%(3*y)], C(0)[i]);
            m += N(0)[i];
            J += X(0)[i+1] * X(1)[i+2] * X(0)[(i+4)%(3*y)] * X(0)[i+2];
            x_field += X(0)[i+1] + X(0)[i+2];
         }
      }
      // TODO
#if 0
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            tx += dot(CH(0)[i], C(0)[i+1]) + dot(C(0)[i], CH(0)[i+1]);
         }
         tx += dot(CH(0)[x-1], C(y+1)[0]) + dot(C(0)[x-1], CH(y+1)[0]);
         for (int i = 0; i < x; ++i)
         {
            ty += dot(CH(0)[i], C(1)[i]) + dot(C(0)[i], CH(1)[i]);
         }
      }
#endif

      Lattice["H_tx"] = sum_unit(tx);
      Lattice["H_ty"] = sum_unit(ty);
      Lattice["H_t"] = sum_unit(tx+ty);
      Lattice["H_m"] = sum_unit(m);
      Lattice["H_J"] = sum_unit(J);
      Lattice["H_x"] = sum_unit(x_field);

      // Gauss's law operators.
      UnitCellMPO G[2*y];

      for (int i = 0; i < y; ++i)
      {
         G[i] = P(0)[3*i] * X(0)[3*i+1] * X(-1)[3*i+1] * X(0)[3*i+2] * X(0)[(3*(i+y)-1)%(3*y)];
         G[i].set_description("Gauss's law operator for site " + std::to_string(3*i));
         Lattice.GetUnitCell().assign_operator("G" + std::to_string(3*i), G[i]);
      }

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
