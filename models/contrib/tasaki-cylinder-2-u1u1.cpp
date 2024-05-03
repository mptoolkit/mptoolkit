// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/tasaki-cylinder-2-u1u1.cpp
//
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1)xU(1) Tasaki model
// This uses a contrived mapping onto an MPS to allow us to generate
// unentangled highest-weight states at quarter filling.
//
// --4--0--1--.--
//      |     |
//      5     .
//      |     |
// --.--3--.--0 (next UC) ...
//      |     |
//      2     .
//      |     |
// --4--0--1--.--
//
//           { 0 p site
// i mod 3 = { 1 horizontal u site
//           { 2 vertical u site

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
      int x = 0;
      int y = 4;
      double nu = 1.0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         //(",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str()) TODO: twisted boundaries
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector (must be even)", y).c_str())
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
      OpDescriptions.set_description("U(1)xU(1) Tasaki cylinder model");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_tx", "nearest-neighbor hopping in y-direction")
         ("H_ty", "nearest-neighbor hopping in x-direction")
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

      if (y % 2 != 0)
      {
         std::cerr << "y must be even" << std::endl;
         return 1;
      }

      int w = x == 0 ? y : x;
      int CellSize = 3*w;

      LatticeSite Site = FermionU1U1();
      UnitCell UCell(Site.GetSymmetryList(), Site, Site, Site);
      UnitCell Cell(repeat(UCell, w));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CHup(Cell, "CHup"), CHdown(Cell, "CHdown"), Cup(Cell, "Cup"),
                       Cdown(Cell, "Cdown"), Pdouble(Cell, "Pdouble"),
                       Hu(Cell, "Hu"), N(Cell, "N"), I(Cell, "I");

      UnitCellMPO AHup[CellSize], AHdown[CellSize], Aup[CellSize], Adown[CellSize];

      for (int i = 0; i < CellSize; i += 3)
      {
         if (i/3 % 2 == 0)
         {
            // Even sites (the A operators will act locally).
            AHup[i] = CHup(0)[i] - nu * (CHup(0)[i+1] + CHup(0)[(i-2+CellSize)%CellSize] + CHup(0)[i+2] + CHup(0)[(i-1+CellSize)%CellSize]) ;
            AHdown[i] = CHdown(0)[i] - nu * (CHdown(0)[i+1] + CHdown(0)[(i-2+CellSize)%CellSize] + CHdown(0)[i+2] + CHdown(0)[(i-1+CellSize)%CellSize]);
            Aup[i] = Cup(0)[i] - nu * (Cup(0)[i+1] + Cup(0)[(i-2+CellSize)%CellSize] + Cup(0)[i+2] + Cup(0)[(i-1+CellSize)%CellSize]);
            Adown[i] = Cdown(0)[i] - nu * (Cdown(0)[i+1] + Cdown(0)[(i-2+CellSize)%CellSize] + Cdown(0)[i+2] + Cdown(0)[(i-1+CellSize)%CellSize]);
         }
         else
         {
         // Odd sites (the A operators will act across the unit cell boundary).
            AHup[i] = CHup(0)[i] - nu * (CHup(1)[(i-5+CellSize)%CellSize] + CHup(-1)[(i+4)%CellSize] + CHup(0)[i+2] + CHup(0)[(i-1+CellSize)%CellSize]) ;
            AHdown[i] = CHdown(0)[i] - nu * (CHdown(1)[(i-5+CellSize)%CellSize] + CHdown(-1)[(i+4)%CellSize] + CHdown(0)[i+2] + CHdown(0)[(i-1+CellSize)%CellSize]);
            Aup[i] = Cup(0)[i] - nu * (Cup(1)[(i-5+CellSize)%CellSize] + Cup(-1)[(i+4)%CellSize] + Cup(0)[i+2] + Cup(0)[(i-1+CellSize)%CellSize]);
            Adown[i] = Cdown(0)[i] - nu * (Cdown(1)[(i-5+CellSize)%CellSize] + Cdown(-1)[(i+4)%CellSize] + Cdown(0)[i+2] + Cdown(0)[(i-1+CellSize)%CellSize]);
         }

         AHup[i].set_description("AHup " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("AHup" + std::to_string(i), AHup[i]);
         AHdown[i].set_description("AHdown " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("AHdown" + std::to_string(i), AHdown[i]);
         Aup[i].set_description("Aup " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("Aup" + std::to_string(i), Aup[i]);
         Adown[i].set_description("Adown " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("Adown" + std::to_string(i), Adown[i]);
      }

      UnitCellMPO BHup[CellSize], BHdown[CellSize], Bup[CellSize], Bdown[CellSize];

      for (int i = 0; i < CellSize; i += 3)
      {
         if (i/3 % 2 == 0)
         {
            // Even sites
            BHup[i+1] = CHup(0)[i+1] + nu * (CHup(0)[i] + CHup(1)[(i-3+CellSize)%CellSize]);
            BHdown[i+1] = CHdown(0)[i+1] + nu * (CHdown(0)[i] + CHdown(1)[(i-3+CellSize)%CellSize]);
            Bup[i+1] = Cup(0)[i+1] + nu * (Cup(0)[i] + Cup(1)[(i-3+CellSize)%CellSize]);
            Bdown[i+1] = Cdown(0)[i+1] + nu * (Cdown(0)[i] + Cdown(1)[(i-3+CellSize)%CellSize]);
         }
         else
         {
            // Odd sites
            BHup[i+1] = CHup(0)[i+1] + nu * (CHup(-1)[(i+6)%CellSize] + CHup(0)[(i+3)%CellSize]);
            BHdown[i+1] = CHdown(0)[i+1] + nu * (CHdown(-1)[(i+6)%CellSize] + CHdown(0)[(i+3)%CellSize]);
            Bup[i+1] = Cup(0)[i+1] + nu * (Cup(-1)[(i+6)%CellSize] + Cup(0)[(i+3)%CellSize]);
            Bdown[i+1] = Cdown(0)[i+1] + nu * (Cdown(-1)[(i+6)%CellSize] + Cdown(0)[(i+3)%CellSize]);
         }

         BHup[i+1].set_description("BHup " + std::to_string(i+1));
         Lattice.GetUnitCell().assign_operator("BHup" + std::to_string(i+1), BHup[i+1]);
         BHdown[i+1].set_description("BHdown " + std::to_string(i+1));
         Lattice.GetUnitCell().assign_operator("BHdown" + std::to_string(i+1), BHdown[i+1]);
         Bup[i+1].set_description("Bup " + std::to_string(i+1));
         Lattice.GetUnitCell().assign_operator("Bup" + std::to_string(i+1), Bup[i+1]);
         Bdown[i+1].set_description("Bdown " + std::to_string(i+1));
         Lattice.GetUnitCell().assign_operator("Bdown" + std::to_string(i+1), Bdown[i+1]);

         BHup[i+2] = CHup(0)[i+2] + nu * (CHup(0)[i] + CHup(0)[(i+3)%CellSize]);
         BHdown[i+2] = CHdown(0)[i+2] + nu * (CHdown(0)[i] + CHdown(0)[(i+3)%CellSize]);
         Bup[i+2] = Cup(0)[i+2] + nu * (Cup(0)[i] + Cup(0)[(i+3)%CellSize]);
         Bdown[i+2] = Cdown(0)[i+2] + nu * (Cdown(0)[i] + Cdown(0)[(i+3)%CellSize]);

         BHup[i+2].set_description("BHup " + std::to_string(i+2));
         Lattice.GetUnitCell().assign_operator("BHup" + std::to_string(i+2), BHup[i+2]);
         BHdown[i+2].set_description("BHdown " + std::to_string(i+2));
         Lattice.GetUnitCell().assign_operator("BHdown" + std::to_string(i+2), BHdown[i+2]);
         Bup[i+2].set_description("Bup " + std::to_string(i+2));
         Lattice.GetUnitCell().assign_operator("Bup" + std::to_string(i+2), Bup[i+2]);
         Bdown[i+2].set_description("Bdown " + std::to_string(i+2));
         Lattice.GetUnitCell().assign_operator("Bdown" + std::to_string(i+2), Bdown[i+2]);
      }

      UnitCellMPO tx, ty, U, Us;
      for (int i = 0; i < CellSize; i += 3)
      {
         tx += dot(BHup[i+1], Bup[i+1]) + dot(BHdown[i+1], Bdown[i+1]);
         ty += dot(BHup[i+2], Bup[i+2]) + dot(BHdown[i+2], Bdown[i+2]);
         U += Pdouble(0)[i] + Pdouble(0)[i+1] + Pdouble(0)[i+2];
         Us += Hu(0)[i] + Hu(0)[i+1] + Hu(0)[i+2];
      }

      Lattice["H_tx"] = sum_unit(tx);
      Lattice["H_ty"] = sum_unit(ty);
      Lattice["H_t"] = sum_unit(tx+ty);
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
