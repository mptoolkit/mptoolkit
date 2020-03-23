// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/toric-code.cpp
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

// Toric code model.
//
// The lattice for a cylinder with w=2 is shown below (bracketed numbers are
// periodically repeated sites).
//
//             +-(3)-+-(7)-+-(11)+
//             |     |     |
//             0     4     8
//             |     |     |
//       +--1--+--5--+--9--+
//       |     |     |     
//       2     6     10
//       |     |     |
// +--3--+--7--+--11-+
//
// Alternatively, rotating by 45 degrees so the unit cell is vertical.
//
// \    /\
// (2)(7) 8
//   \/    \/
//   /\    /\
// (3) 4  9
// /    \/
// \    /\
//  0  5  10
//   \/    \/
//   /\    /\
//  1  6  11
// /    \/
// \    /\
//  2  7 (8)
//   \/    \/
//   /\    /\
//  3 (4)(9)
// /    \/

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin.h"
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
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the spin", Spin).c_str())
         ("width,w", prog_opt::value(&w), FormatDefault("width of the cylinder", w).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ("noreflect", prog_opt::bool_switch(&NoReflect),
          "don't include the spatial reflection operator (expensive for large width lattices)")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.description("Toric code");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("WX"         , "Wilson loop of X operators across the unit cell")
         ("WZ"         , "Wilson loop of Z operators across the unit cell")
         ("Trans"      , "translation by one site (rotation by 2\u0071/w) in lattice short direction")
         ("Ref"        , "reflection in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("RyUnit"     , "reflection of a single unit cell",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ;
      OpDescriptions.add_operators()
         ("H_x"        , "magnetic field in x direction")
         ("H_y"        , "magnetic field in y direction")
         ("H_z"        , "magnetic field in z direction")
         ("H_star"     , "sum of star operators")
         ("H_plaq"     , "sum of plaquette operators")
         ("Ty"         , "momentum operator in lattice short direction")
         ("TyPi"       , "translation by w sites in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("Ry"         , "reflection in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
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

      LatticeSite Site = SpinSite(Spin);
      UnitCell Cell(repeat(Site, u));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator I(Cell, "I"), X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");
      UnitCellOperator WX(Cell, "WX"), WZ(Cell, "WZ");
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");
      
      // Magnetic fields.
      UnitCellMPO H_x, H_y, H_z;

      for (int i = 0; i < u; ++i)
      {
         H_x += X(0)[i];
         H_y += Y(0)[i];
         H_z += Z(0)[i];
      }

      Lattice["H_x"] = sum_unit(H_x);
      Lattice["H_y"] = sum_unit(H_y);
      Lattice["H_z"] = sum_unit(H_z);

      // Star and plaquette operators.
      UnitCellMPO A[w], B[w];

      for (int i = 0; i < w; ++i)
      {
         A[i] = X(0)[2*i] * X(1)[(2*i+1)%u] * X(0)[(2*i+1)%u] * X(1)[(2*i+2)%u];
         B[i] = Z(0)[(2*i+1)%u] * Z(1)[(2*i+2)%u] * Z(0)[(2*i+2)%u] * Z(1)[(2*i+3)%u];
         A[i].set_description("star operator " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("A" + std::to_string(i), A[i]);
         B[i].set_description("plaquette operator " + std::to_string(i));
         Lattice.GetUnitCell().assign_operator("B" + std::to_string(i), B[i]);
      }

      // Sum of the star and plaquette operators.
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
      for (int i = 0; i < u; i++)
      {
         WX = WX * X(0)[i];
         WZ = WZ * Z(0)[i];
      }

      // Translation and relfection operators.
      Trans = I(0);
      for (int i = 0; i < u-1; ++i)
      {
         //T *= 0.5*( 0.25*inner(S[i],S[i+1]) + 1 );
         Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
      }

      if (!NoReflect)
      {
         Ref = I(0); // old way of representing an explicit R-operator.
         for (int i = 0; i < w; ++i)
         {
            //R *= 0.5*( 0.25*inner(S[i],S[w-i-1]) + 1 );
            Ref = Ref(0) * Cell.swap_gate_no_sign(i, u-i-1);
         }
      }

      UnitCellMPO Ry = I(0);
      if (!NoReflect)
      {
         for (int c = 0; c < u; ++c)
         {
            UnitCellMPO ThisR = I(0);
            // get the 'pivot' site/bond that we reflect about
            int const p1 = c/2;
            int const p2 = (c+1)/2;

            // if we're reflecting about a bond, do that first
            if (p1 != p2)
               ThisR = ThisR * Cell.swap_gate_no_sign(p1,p2);

            int i1 = (p1+u-1)%u;
            int i2 = (p2+1)%u;

            while (i1 != p1 + w)
            {
               ThisR = ThisR * Cell.swap_gate_no_sign(i1,i2);
               i1 = (i1+u-1)%u;
               i2 = (i2+1)%u;
            }

            ThisR.translate(c*u);
            Ry = Ry * ThisR;
         }
         RyUnit = Ry;
      }

      // Momentum operators in Y-direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), u);

      if (!NoReflect)
         Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), u*u);

      // add rotation by pi
      if (!NoReflect)
      {
         UnitCellMPO TyPi = I(0);
         for (int i = 0; i < w; ++i)
         {
            TyPi = TyPi * Cell.swap_gate_no_sign(i, i+w);
         }
         Lattice["TyPi"] = prod_unit_left_to_right(TyPi.MPO(), u);
      }

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
