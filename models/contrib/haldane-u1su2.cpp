// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/haldane-u1su2.cpp
//
// Copyright (C) 2020 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1)xSU(2) Haldane model.
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
// Nearest-neighbour interactions are:
// (vertical)    (0)[i] (0)[(i+1)%(2w)]  (for i = 0..2w), 2w terms per unit cell
// (horizontal)  (0)[i] (1)[(i+1)%(2w)]  (for i = 0,2,4,...,2(w-1)), w terms per unit cell

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/fermion-u1su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

using math_const::pi;
namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int w = 4;
      std::string FileName;
      bool NoReflect = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
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
      OpDescriptions.set_description("U(1)xSU(2) Haldane model");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("Trans"      , "translation by one site (rotation by 2\u0071/w) in lattice short direction")
         ("Ref"        , "reflection in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("RyUnit"     , "reflection of a single unit cell",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ;
      OpDescriptions.add_operators()
         ("H_t1"       , "nearest-neighbour hopping")
         ("H_t2cw"     , "next-nearest-neighbour clockwise hopping")
         ("H_t2acw"    , "next-nearest-neighbour anticlockwise hopping")
         ("H_t3"       , "third-nearest-neighbour hopping")
         ("H_M"        , "asymmetric potential energy")
         ("H_U"        , "on-site interaction n_up*n_down")
         ("H_Us"       , "on-site interaction (n_up-1/2)*(n_down-1/2)")
         ("H_V1"       , "nearest-neighbour interaction")
         ("H_V2"       , "next-nearest-neighbour interaction")
         ("Ty"         , "momentum operator in lattice short direction")
         ("TyPi"       , "translation by w sites in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("Ry"         , "reflection in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ;
      OpDescriptions.add_functions()
         ("H_t2"       , "next-nearest-neighbour complex hopping with phase phi")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = FermionU1SU2();
      int u = w*2;
      UnitCell Cell = repeat(Site, u);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator I(Cell, "I"), CH(Cell, "CH"), C(Cell, "C"), Pdouble(Cell, "Pdouble"), Hu(Cell, "Hu"), N(Cell, "N"), N2(Cell, "N2");
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      // Haldane model terms.
      UnitCellMPO H_t1, H_t2cw, H_t2acw, H_t3, H_M;

      for (int i = 0; i < u; i += 2)
      {
         H_t1 += dot(CH(0)[i], C(0)[(i+1)%u]) + dot(C(0)[i], CH(0)[(i+1)%u])
               + dot(CH(0)[i], C(1)[(i+1)%u]) + dot(C(0)[i], CH(1)[(i+1)%u])
               + dot(CH(0)[i], C(0)[(i+u-1)%u]) + dot(C(0)[i], CH(0)[(i+u-1)%u]);
         H_t2cw += dot(CH(0)[i], C(0)[(i+2)%u])
                 + dot(CH(0)[(i+1)%u], C(1)[(i+3)%u])
                 + dot(CH(0)[(i+2)%u], C(1)[(i+2)%u])
                 + dot(CH(1)[(i+3)%u], C(1)[(i+1)%u])
                 + dot(CH(1)[(i+2)%u], C(0)[i])
                 + dot(CH(1)[(i+1)%u], C(0)[(i+1)%u]);
         H_t2acw += dot(C(0)[i], CH(0)[(i+2)%u])
                  + dot(C(0)[(i+1)%u], CH(1)[(i+3)%u])
                  + dot(C(0)[(i+2)%u], CH(1)[(i+2)%u])
                  + dot(C(1)[(i+3)%u], CH(1)[(i+1)%u])
                  + dot(C(1)[(i+2)%u], CH(0)[i])
                  + dot(C(1)[(i+1)%u], CH(0)[(i+1)%u]);
         H_t3 += dot(CH(0)[i], C(1)[(i+3)%u]) + dot(C(0)[i], CH(1)[(i+3)%u])
               + dot(CH(0)[i], C(1)[(i+u-1)%u]) + dot(C(0)[i], CH(1)[(i+u-1)%u])
               + dot(CH(0)[i], C(-1)[(i+u-1)%u]) + dot(C(0)[i], CH(-1)[(i+u-1)%u]);
         H_M += N(0)[i] - N(0)[(i+1)%u];
      }

      Lattice["H_t1"] = sum_unit(H_t1);
      Lattice["H_t2cw"] = sum_unit(H_t2cw);
      Lattice["H_t2acw"] = sum_unit(H_t2acw);
      Lattice.func("H_t2")(arg("phi")) = "exp(i*phi)*H_t2cw + exp(-i*phi)*H_t2acw";
      Lattice["H_t3"] = sum_unit(H_t3);
      Lattice["H_M"] = sum_unit(H_M);

      // On-site interaction.
      UnitCellMPO H_U, H_Us;

      for (int i = 0; i < u; ++i) {
         H_U += Pdouble(0)[i];
         H_Us += Hu(0)[i];
      }
      
      Lattice["H_U"] = sum_unit(H_U);
      Lattice["H_Us"] = sum_unit(H_Us);

      // Neighbouring site interaction.
      UnitCellMPO H_V1, H_V2;

      for (int i = 0; i < u; i += 2)
      {
         H_V1 += dot(N(0)[i], N(0)[(i+1)%u])
               + dot(N(0)[i], N(1)[(i+1)%u])
               + dot(N(0)[i], N(0)[(i+u-1)%u]);
         H_V2 += dot(N(0)[i], N(0)[(i+2)%u])
               + dot(N(0)[(i+1)%u], N(1)[(i+3)%u])
               + dot(N(0)[(i+2)%u], N(1)[(i+2)%u])
               + dot(N(1)[(i+3)%u], N(1)[(i+1)%u])
               + dot(N(1)[(i+2)%u], N(0)[i])
               + dot(N(1)[(i+1)%u], N(0)[(i+1)%u]);
      }

      Lattice["H_V1"] = sum_unit(H_V1);
      Lattice["H_V2"] = sum_unit(H_V2);

      // Translation and reflection operators.
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
