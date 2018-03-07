// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-hex-zigzag-su2.cpp
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// Hexagonal lattice spin model
// Example for a width-6 lattice (site numbers in brackets are periodic repeats
// in the vertical
// direction (i.e. top-left (5) is the same site as the bottom-left 5).
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

// Interactions are:
// (vertical)    (0)[i] (0)[(i+1)%w]  (for i = 0..w), w terms per unit cell
// (horizontal)  (0)[i] (1)[(i+1)%w]  (for i = 0,2,4,...,w-1), (w/2) terms per unit cell

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

using math_const::pi;
namespace prog_opt = boost::program_options;


std::complex<double> phase(double theta)
{
   return std::exp(std::complex<double>(0.0, theta));
}


int IntPow(int x, int p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * IntPow(x, p-1);
}

// Functions for producing dimer terms in the Hamiltonian (including their h.c. (rotated) terms):
UnitCellMPO DimerKinetic(half_int spin, int width, int cell1, int site1,
                         int cell2, int site2, int cell3, int site3, int cell4, int site4) {

 LatticeSite Site = SpinSU2(spin);
 UnitCell Cell = repeat(Site, width);
 UnitCellOperator S(Cell, "S"), I(Cell, "I");
 UnitCellMPO H, Dots, DoubleDots, CrossDots;

 Dots = inner(S(cell1)[site1], S(cell4)[site4]) + inner(S(cell2)[site2], S(cell3)[site3])
    - inner(S(cell1)[site1], S(cell2)[site2])
    - inner(S(cell3)[site3], S(cell4)[site4]) - inner(S(cell1)[site1], S(cell3)[site3])
    - inner(S(cell2)[site2], S(cell4)[site4]);

 DoubleDots = inner(S(cell1)[site1], S(cell2)[site2])*inner(S(cell3)[site3], S(cell4)[site4]) +
              inner(S(cell1)[site1], S(cell3)[site3])*inner(S(cell2)[site2], S(cell4)[site4]);

 CrossDots = inner(cross(S(cell1)[site1], S(cell2)[site2]), cross(S(cell3)[site3], S(cell4)[site4])) +
             inner(cross(S(cell1)[site1], S(cell3)[site3]), cross(S(cell2)[site2], S(cell4)[site4]));

 H = 0.125*I + 0.25*Dots + 0.5*DoubleDots - 0.125*CrossDots;

 return H;

}


int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 3;
      std::string FileName;
      bool NoReflect = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("width,w", prog_opt::value(&w), "width of the cylinder [default 3]")
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
      OpDescriptions.set_description("SU(2) hexagonal lattice zigzag configuration");
      OpDescriptions.author("IP McCullocch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_cell_operators()
         ("S"          , "total spin on a leg of the cylinder")
         ("Trans"      , "translation by one site (rotation by 2\u0071/w) in lattice short direction")
         ("Ref"        , "reflection in lattice short direction",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ("RyUnit"     , "Reflection of a single unit cell",
          "not present with --noreflect", [&NoReflect]()->bool{return !NoReflect;})
         ;
      OpDescriptions.add_operators()
         ("H_J1"       , "nearest neighbor spin exchange")
         ("H_y"        , "nearest neighbor spin exchange in the Y direction")
         ("Ty"         , "momentum operator in lattice short direction")
         ("TyPi"       , "translation by w/2 sites in the Y direction",
          "width even, not present with --noreflect",
          [&NoReflect,&w]()->bool{return !NoReflect && w%2 == 0;})
         ("Ry"         , "Reflection in the Y direction",
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

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, w);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator S(Cell, "S");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      for (int i = 0; i < w; ++i)
      {
         S += S[i];                    // total spin on a leg of cylinder
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
       {
           //T *= 0.5*( 0.25*inner(S[i],S[i+1]) + 1 );
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
       }

      if (!NoReflect)
      {
         Ref = I(0); // old way of representing an explicit R-operator.
         for (int i = 0; i < w/2; ++i)
         {
            //R *= 0.5*( 0.25*inner(S[i],S[w-i-1]) + 1 );
            Ref = Ref(0) * Cell.swap_gate_no_sign(i, w-i-1);
         }
      }


      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2, Ht, Hv, Hy;

      for (int i = 0; i < w; ++i)
      {
        // nearest neighbor bonds
         H1 += inner(S(0)[i], S(0)[(i+1)%w]);   // vertical bonds
	 if (i % 2 == 0)
	    H1 += inner(S(0)[i], S(1)[(i+1)%w]);

         // vertical bonds only
         Hy += inner(S(0)[i], S(0)[(i+1)%w]);

      }

      // Reflection.  This is in the 'wrong' 45 degree angle
      UnitCellMPO Ry = I(0);
      if (!NoReflect)
      {
         for (int c = 0; c < w; ++c)
         {
            UnitCellMPO ThisR = I(0);
            // get the 'pivot' site/bond that we reflect about
            int const p1 = c/2;
            int const p2 = (c+1)/2;

            // if we're reflecting about a bond, do that first
            if (p1 != p2)
               ThisR = ThisR * Cell.swap_gate_no_sign(p1,p2);

            int i1 = (p1+w-1)%w;
            int i2 = (p2+1)%w;

            while (i1 != p1 + w/2)
            {
               ThisR = ThisR * Cell.swap_gate_no_sign(i1,i2);
               i1 = (i1+w-1)%w;
               i2 = (i2+1)%w;
            }

            ThisR.translate(c*w);
            Ry = Ry * ThisR;
         }
         RyUnit = Ry;
      }

      // Now we construct the InfiniteLattice,

      Lattice["H_J1"]  = sum_unit(H1);
      Lattice["H_y"]   = sum_unit(Hy);

      // Momentum operators in Y-direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);

      if (!NoReflect)
         Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w*w);

      // for even size unit cell, add rotation by pi
      if (w%2 == 0 && !NoReflect)
      {
         UnitCellMPO TyPi = I(0);
         for (int i = 0; i < w/2; ++i)
         {
            TyPi = TyPi * Cell.swap_gate_no_sign(i, i+w/2);
         }
         Lattice["TyPi"] = prod_unit_left_to_right(TyPi.MPO(), w);
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
