// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/u1-qlm-cylinder-hex-u1.cpp
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

// U(1) 2+1D U(1) quantum link model on a hexagonal (honeycomb) lattice
//
//          0
//         /
//        8
//       /
// --9--6            (5)
//       \           /
//        7        (3)
//         \       /
//          5-(4)(1)
//         /       \
//        3        (2)
//       /           \
// --4--1            (0)
//       \           /
//        2        (8)
//         \       /
//          0-(9)(6)
//
//      Y
//     /
// Z--+
//     \
//      X
//
// 0 and 5 are fermions, while 1 and 6 are antifermions.
// The other sites are link sites represented by spins.

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
      int y = 4;
      half_int Spin = 0.5;
      bool Bosonic = false;
      bool OBC = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector (must be even)", y).c_str())
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the link spins", Spin).c_str())
         ("bosonic", prog_opt::bool_switch(&Bosonic), "use hardcore bosons instead of spinless fermions on the matter sites")
         ("obc", prog_opt::bool_switch(&OBC), "use open boundary conditions along the y direction")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) 2+1D U(1) quantum link model hexagonal lattice (U.-J. Wiese, Annalen der Physik 525, 777 (2013), arXiv:1305.1602)");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_tx" , "nearest-neighbor hopping in x-direction")
         ("H_ty" , "nearest-neighbor hopping in y-direction")
         ("H_tz" , "nearest-neighbor hopping in z-direction")
         ("H_t"  , "nearest-neighbor hopping")
         ("H_m"  , "fermion mass")
         ("H_g"  , "gauge coupling")
         ("H_chi"      , "background field")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         std::cerr << "Operators:" << std::endl << OpDescriptions;
         return 1;
      }

      if (y % 2 != 0)
      {
         std::cerr << "y must be even" << std::endl;
         return 1;
      }

      int CellSize = 5 * (y/2);

      LatticeSite FSite = SpinlessFermionU1("N", "P", Bosonic);
      LatticeSite AFSite = SpinlessAntifermionU1("N", "P", Bosonic);
      LatticeSite SSite = SpinSite(Spin);

      UnitCell NCell(FSite.GetSymmetryList(), FSite, AFSite);
      UnitCell SCell(FSite.GetSymmetryList(), SSite, SSite, SSite);
      UnitCell Cell(repeat(join(NCell, SCell), CellSize/5));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N"), I(Cell, "I"),
                       Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");

      UnitCellMPO tx, ty, tz, m, g, chi;

      for (int i = 0; i < CellSize; i += 5)
      {
         tx += Sp(0)[i+2] * dot(CH(0)[i],   C(0)[i+1])             + Sm(0)[i+2] * dot(CH(0)[i+1],             C(0)[i]);
         if (!(OBC && i == CellSize - 5)) // Do not include these terms for the final site if using OBC
         {
            ty += Sp(0)[i+3] * dot(CH(0)[i+1], C(0)[(i+5)%CellSize])  + Sm(0)[i+3] * dot(CH(0)[(i+5)%CellSize],  C(0)[i+1]);
            tz += Sp(0)[i+4] * dot(CH(0)[i+1], C(-1)[(i+5)%CellSize]) + Sm(0)[i+4] * dot(CH(-1)[(i+5)%CellSize], C(0)[i+1]);
         }
         chi += Sz(0)[i+2] + Sz(0)[i+3] + Sz(0)[i+4];
         g += Sz(0)[i+2]*Sz(0)[i+2] + Sz(0)[i+3]*Sz(0)[i+3] + Sz(0)[i+4]*Sz(0)[i+4];
         m += N(0)[i] - N(0)[i+1];
      }

      Lattice["H_tx"] = sum_unit(tx);
      Lattice["H_ty"] = sum_unit(ty);
      Lattice["H_tz"] = sum_unit(tz);
      Lattice["H_t"] = sum_unit(tx+ty+tz);
      Lattice["H_m"] = sum_unit(m);
      Lattice["H_g"] = sum_unit(g);
      Lattice["H_chi"] = sum_unit(chi);

      // Gauss's law operators.
      UnitCellMPO G[y];

      for (int j = 0; j < y; j += 2)
      {
         int i = 5*(j/2);

         G[j] = N(0)[i] - Sz(0)[i+2] + Sz(0)[(i-2+CellSize)%CellSize] + Sz(1)[(i-1+CellSize)%CellSize];
         G[j].set_description("Gauss's law operator for site " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("G" + std::to_string(i), G[j]);

         G[j+1] = N(0)[i+1] + Sz(0)[i+2] - Sz(0)[i+3] - Sz(0)[i+4];
         G[j+1].set_description("Gauss's law operator for site " + std::to_string(i+1));
         Lattice.GetUnitCell().assign_operator("G" + std::to_string(i+1), G[j+1]);
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
