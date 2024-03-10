// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/toric-code-z2.cpp
//
// Copyright (C) 2020 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// Z2 Toric code model.
//
// The lattice for a cylinder with w=2 is shown below (the bracketed numbers
// represent periodically repeated sites).
//
// +-(3)-+-(7)-+-(11)
// |     |     |
// 0     4     8
// |     |     |
// +--1--+--5--+--9--
// |     |     |
// 2     6     10
// |     |     |
// +--3--+--7--+--11-

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-z2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 4;
      std::string FileName;
      bool NoReflect = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         //("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the spin", Spin).c_str()) // This should only be 0.5.
         ("width,w", prog_opt::value(&w), FormatDefault("width of the cylinder", w).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.description("Z2 Toric code");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("WX"         , "Wilson loop of X operators around the circumference")
         ("WZ"         , "Wilson loop of Z operators around the circumference")
         ("Trans"      , "translation by one site (rotation by 2\u0071/w) in lattice short direction")
         ("Ref"        , "reflection in lattice short direction")
         ;
      OpDescriptions.add_operators()
         ("H_x"        , "magnetic field in x direction")
         ("H_star"     , "sum of star operators")
         ("H_plaq"     , "sum of plaquette operators")
         ("Ty"         , "momentum operator in lattice short direction")
         ("TyPi"       , "translation by w sites in lattice short direction")
         ("Ry"         , "reflection in lattice short direction")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cout << OpDescriptions << '\n';
         return 1;
      }

      int u = 2*w;

      LatticeSite Site = SpinZ2(Spin);
      UnitCell Cell(repeat(Site, u));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator I(Cell, "I"), X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");
      UnitCellOperator WX(Cell, "WX"), WZ(Cell, "WZ");
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");
      
      // Magnetic field.
      UnitCellMPO H_x;

      for (int i = 0; i < u; ++i)
         H_x += X(0)[i];

      Lattice["H_x"] = sum_unit(H_x);

      // Star and plaquette operators.
      UnitCellMPO A[w], B[w];

      for (int i = 0; i < w; ++i)
      {
         A[i] = X(1)[2*i] * X(1)[2*i+1] * X(0)[2*i+1] * X(1)[(2*i+2)%u];
         B[i] = Z(0)[(2*i+u-1)%u] * Z(0)[2*i] * Z(1)[2*i] * Z(0)[2*i+1];
         A[i].set_description("star operator " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("A" + std::to_string(i), A[i]);
         B[i].set_description("plaquette operator " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("B" + std::to_string(i), B[i]);
      }

      // Sums of the star and plaquette operators.
      UnitCellMPO H_star, H_plaq;

      for (int i = 0; i < w; ++i)
      {
         H_star += A[i];
         H_plaq += B[i];
      }

      Lattice["H_star"] = sum_unit(H_star);
      Lattice["H_plaq"] = sum_unit(H_plaq);

      // Wilson loop operators.
      WX = I(0);
      WZ = I(0);
      for (int i = 0; i < w; ++i)
      {
         WX = WX * X(0)[2*i+1];
         WZ = WZ * Z(0)[2*i];
      }

      // Translation operators.
      Trans = I(0);
      for (int i = 0; i < u-2; ++i)
         Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+2);

      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), u);

      // Reflection operators.
      Ref = I(0);
      for (int i = 1; i < w; ++i)
         Ref = Ref(0) * Cell.swap_gate_no_sign(i, u-i);

      Lattice["Ry"] = prod_unit_left_to_right(UnitCellMPO(Ref(0)).MPO(), u);

      // Rotation by pi.
      UnitCellMPO TyPi = I(0);
      for (int i = 0; i < w; ++i)
         TyPi = TyPi * Cell.swap_gate_no_sign(i, i+w);

      Lattice["TyPi"] = prod_unit_left_to_right(TyPi.MPO(), u);

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
