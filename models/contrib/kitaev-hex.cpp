// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/kitaev-hex.cpp
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

// Kitaev honeycomb model
// Horizontal bonds are Sx, +60 degree bonds are Sy, -60 degree bonds are Sz.
//
// Example for a width-3 lattice (site numbers in brackets are periodic repeats
// in the vertical direction (i.e. top-left (5) is the same site as the
// bottom-left 5).
// Sites 6,7,8,9,10,11 are the second unit cell, e.g. 6 is (1)[0].
//
// (4)-(11)    12
//  /    \    /
//(5)     6--13
//  \    /    \
//   0--7      14
//  /    \    /
// 1      8--15
//  \    /    \
//   2--9      16
//  /    \    /
// 3      10-17
//  \    /    \
//   4--11    (12)
//  /    \    /
// 5     (6)-(13)
//  \    /
//  (0)-(7)
//
//     Y
//    /
// X--
//    \
//     Z
// Interactions are:
// (vertical)    (0)[i] (0)[(i+1)%(2w)]  (for i = 0..2w), 2w terms per unit cell
// (horizontal)  (0)[i] (1)[(i+1)%(2w)]  (for i = 0,2,4,...,2(w-1)), w terms per unit cell

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

using math_const::pi;
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
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("width,w", prog_opt::value(&w), "width of the cylinder [default 4]")
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

      // Descriptions of each operator
      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("Kitaev honeycomb model");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("WY"         , "Wilson loop in the Y direction")
         ("Trans"      , "translation by one site (rotation by 2\u0071/w) in lattice short direction")
         ("Ref"        , "reflection in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("RyUnit"     , "reflection of a single unit cell",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ;
      OpDescriptions.add_operators()
         ("H_x"        , "magnetic field in the x direction")
         ("H_y"        , "magnetic field in the y direction")
         ("H_z"        , "magnetic field in the z direction")
         ("H_xx"       , "horizontal X-X interaction")
         ("H_yy"       , "+60 degrees Y-Y interaction")
         ("H_zz"       , "-60 degrees Z-Z interaction")
         ("H_3"        , "three-spin interaction")
         ("Ty"         , "momentum operator in lattice short direction")
         ("TyPi"       , "translation by w sites in the Y direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("Ry"         , "reflection in the Y direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = SpinSite(Spin);
      int u = w*2;
      UnitCell Cell = repeat(Site, u);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator I(Cell, "I"), X(Cell, "X"), Y(Cell, "Y"), Z(Cell, "Z");
      UnitCellOperator WY(Cell, "WY");
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      // Plaquette operators.
      for (int i = 0; i < u; i += 2)
      {
         UnitCellMPO W = Z(0)[i] * X(0)[(i+1)%u] * Y(0)[(i+2)%u]
            * Z(1)[(i+3)%u] * X(1)[(i+2)%u] * Y(1)[(i+1)%u];
         W.set_description("plaquette operator " + std::to_string(i/2));
         std::string Name = "W" + std::to_string(i/2);
         Lattice.GetUnitCell().assign_operator(Name, W);
      }

      // Wilson loop operator.
      WY = I(0);

      for (int i = 0; i < u; ++i)
      {
         WY = WY * X(0)[i];
      }

      // Magnetic fields.
      UnitCellMPO Hx, Hy, Hz;

      for (int i = 0; i < u; ++i)
      {
         Hx += X(0)[i];
         Hy += Y(0)[i];
         Hz += Z(0)[i];
      }

      Lattice["H_x"] = sum_unit(Hx);
      Lattice["H_y"] = sum_unit(Hy);
      Lattice["H_z"] = sum_unit(Hz);

      // Kitaev model interactions.
      UnitCellMPO Hxx, Hyy, Hzz;

      for (int i = 0; i < u; i += 2)
      {
         Hxx -= X(0)[i] * X(1)[(i+1)%u];
         Hyy -= Y(0)[i] * Y(0)[(i+1)%u];
         Hzz -= Z(0)[i] * Z(0)[(i+u-1)%u];
      }

      Lattice["H_xx"] = sum_unit(Hxx);
      Lattice["H_yy"] = sum_unit(Hyy);
      Lattice["H_zz"] = sum_unit(Hzz);

      // Three spin interaction term.
      UnitCellMPO H3;

      for (int i = 0; i < u; i += 2)
      {
         H3 += Y(0)[(i+1)%u] * X(1)[(i+1)%u] * Z(0)[i]
             + X(0)[(i+1)%u] * Z(0)[(i+2)%u] * Y(0)[i]
             + Z(0)[(i+1)%u] * Y(0)[(i+2)%u] * X(1)[(i+3)%u]
             + Y(1)[(i+2)%u] * X(0)[(i+2)%u] * Z(1)[(i+3)%u]
             + X(1)[(i+2)%u] * Z(1)[(i+1)%u] * Y(1)[(i+3)%u]
             + Z(1)[(i+2)%u] * Y(1)[(i+1)%u] * X(0)[i];
      }

      Lattice["H_3"] = sum_unit(H3);

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

      // Reflection.
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
