// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/u1-qlm-cylinder-u1.cpp
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1) 2+1D U(1) quantum link model
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
      int x = 0;
      int y = 4;
      half_int Spin = 0.5;
      bool Bosonic = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         //(",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str()) TODO: twisted boundaries
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector (must be even)", y).c_str())
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the link spins", Spin).c_str())
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
      OpDescriptions.set_description("U(1) 2+1D U(1) quantum link model (U.-J. Wiese, Annalen der Physik 525, 777 (2013), arXiv:1305.1602)");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_tx" , "nearest-neighbor hopping in y-direction")
         ("H_ty" , "nearest-neighbor hopping in x-direction")
         ("H_t"  , "nearest-neighbor hopping")
         ("H_m"  , "fermion mass")
         ("H_g"  , "gauge coupling")
         ("H_J"  , "plaquette interactions")
         ("H_flux"     , "electric flux")
         ("H_stag_flux", "staggered electric flux")
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

      int CellSize = x == 0 ? 6*y : 6*x;

      LatticeSite FSite = SpinlessFermionU1("N", "P", Bosonic);
      LatticeSite AFSite = SpinlessAntifermionU1("N", "P", Bosonic);
      LatticeSite SSite = SpinSite(Spin);
      UnitCell FCell(FSite.GetSymmetryList(), FSite, SSite, SSite);
      UnitCell AFCell(AFSite.GetSymmetryList(), AFSite, SSite, SSite);
      UnitCell Cell1(repeat(join(FCell, AFCell), CellSize/12));
      UnitCell Cell2(repeat(join(AFCell, FCell), CellSize/12));
      UnitCell Cell = join(Cell1, Cell2);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N"), I(Cell, "I"),
                       Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");

      UnitCellMPO tx, ty, m, g, J, flux, stag_flux;
      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < 3*y; i += 3)
         {
            tx += Sp(0)[i+1]     * dot(CH(0)[i],     C(0)[i+3*y])           + Sm(0)[i+1]     * dot(CH(0)[i+3*y], C(0)[i]);
            tx += Sp(0)[i+1+3*y] * dot(CH(0)[i+3*y], C(1)[i])               + Sm(0)[i+1+3*y] * dot(CH(1)[i], C(0)[i+3*y]);
            ty += Sp(0)[i+2]     * dot(CH(0)[i],     C(0)[(i+3)%(3*y)])     + Sm(0)[i+2]     * dot(CH(0)[(i+3)%(3*y)], C(0)[i]);
            // This term is intentionally negative.
            ty -= Sp(0)[i+2+3*y] * dot(CH(0)[i+3*y], C(0)[(i+3)%(3*y)+3*y]) + Sm(0)[i+2+3*y] * dot(CH(0)[(i+3)%(3*y)+3*y], C(0)[i+3*y]);
            J += Sp(0)[i+1] * Sp(0)[i+2+3*y] * Sm(0)[(i+4)%(3*y)] * Sm(0)[i+2];
            J += Sm(0)[i+1] * Sm(0)[i+2+3*y] * Sp(0)[(i+4)%(3*y)] * Sp(0)[i+2];
            J += Sp(0)[i+1+3*y] * Sp(1)[i+2] * Sm(0)[(i+4)%(3*y)+3*y] * Sm(0)[i+2+3*y];
            J += Sm(0)[i+1+3*y] * Sm(1)[i+2] * Sp(0)[(i+4)%(3*y)+3*y] * Sp(0)[i+2+3*y];
            flux += Sz(0)[i+1] + Sz(0)[i+2] + Sz(0)[i+1+3*y] + Sz(0)[i+2+3*y];
            g += Sz(0)[i+1]*Sz(0)[i+1] + Sz(0)[i+2]*Sz(0)[i+2] + Sz(0)[i+1+3*y]*Sz(0)[i+1+3*y] + Sz(0)[i+2+3*y]*Sz(0)[i+2+3*y];
         }
         for (int i = 0; i < 3*y; i += 6)
         {
            m += N(0)[i] - N(0)[i+3] - N(0)[i+3*y] + N(0)[i+3+3*y];
            stag_flux += 4.0*I(0)[i+1] - Sz(0)[i+1] - Sz(0)[i+2] + Sz(0)[i+4] + Sz(0)[i+5] + Sz(0)[i+1+3*y] + Sz(0)[i+2+3*y] - Sz(0)[i+4+3*y] - Sz(0)[i+5+3*y];
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
      Lattice["H_g"] = sum_unit(g);
      Lattice["H_J"] = sum_unit(J);
      Lattice["H_flux"] = sum_unit(flux);
      Lattice["H_stag_flux"] = sum_unit(stag_flux);

      // Gauss's law operators.
      UnitCellMPO G[2*y];

      for (int i = 0; i < y; ++i)
      {
         G[i] = N(0)[3*i] - Sz(0)[3*i+1] + Sz(-1)[3*(i+y)+1] - Sz(0)[3*i+2] + Sz(0)[(3*(i+y)-1)%(3*y)];
         G[i].set_description("Gauss's law operator for site " + std::to_string(3*i));
         Lattice.GetUnitCell().assign_operator("G" + std::to_string(3*i), G[i]);
      }

      for (int i = y; i < 2*y; ++i)
      {
         G[i] = N(0)[3*i] - Sz(0)[3*i+1] + Sz(0)[3*(i-y)+1] - Sz(0)[3*i+2] + Sz(0)[(3*i-1)%(3*y)+3*y];
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
