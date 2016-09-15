// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-tri-yc-z2.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2015-2016 Seyed N. Saadatmand <s.saadatmand@uq.edu.au>
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

// YC configuration of a triangular lattice.
// The default unit-cell size is the width value.
//
// Example for a width-6 lattice (site numbers in brackets are periodic repeats in the vertical
// direction (i.e. top-left (5) is the same site as the bottom-left 5).
// Sites 6,7,8,9,10,11 are the second unit cell, e.g. 6 is (1)[0].
//
//      (17)
//  (11) |
//(5) | 12
// |  6< |
// 0< | 13
// |  7< |
// 1< | 14
// |  8< | ...
// 2< | 15
// |  9< |
// 3< | 16
// | 10< |
// 4< | 17
// | 11< |
// 5< | (12)
// | (6)
//(0)


#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-z2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>


namespace prog_opt = boost::program_options;


int IntPow(int x, int p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * IntPow(x, p-1);
}


std::ostream &OpProgress(int *n, int step, int max)
{
    (*n) = (*n) + step;
    std::cout << "\33[2K\r" << "working... %" << (100*(*n))/max << std::flush;
    return std::cout;
}


int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 3;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("width,w", prog_opt::value(&w), "width of the cylinder [default 3]")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Constructs a triangular lattice in the YC configuration with wrapping vector (0,1),\n"
                   << "while spins sit on edges. This employs an efficient way of numbering in the 1D chain.\n\n"
                   << "Operators:\n"
                   << "Sx      - total spin on a leg of the cylinder\n"
                   << "Sy      - total spin on a leg of the cylinder\n"
                   << "Sz      - total spin on a leg of the cylinder\n"
                   << "Trans   - translation by one site (rotation by 2\u0071/w) in lattice short direction\n"
                   << "Ref     - reflection in lattice short direction (may need applying T-operators to become gneral reflection)\n"
                   << "Functions:\n"
                   << "TIM{J = NN coupling strength, Gamma = strength of magnetic field term in X-direcion}\n\n"
            ;
         return 1;
      }

      const double PI = 4.0*std::atan(1.0);
      //const std::complex<double> jj(0.0,theta*PI);
      int oo = 0;
      int oo_max = 8*w+(w*w/2)+13;

      LatticeSite Site = SpinZ2(Spin);
      UnitCell Cell = repeat(Site, w);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      std::cout << "Building all Hamiltonian operators:\n";

      // Add some operators on the unit-cell

      for (int i = 0; i < w; ++i)
      {
         Sx += Sx[i];                                         // total spin (x-component) on a leg of cylinder
         Sy += Sy[i];                                         // total spin (y-component) on a leg of cylinder
         Sz += Sz[i];                                         // total spin (z-component) on a leg of cylinder

         //oo+=3;
         //std::printf("\33[2K\r");
         //std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 3*w
         OpProgress(&oo,3,oo_max); // operator series count: 3*w
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
       {
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
            OpProgress(&oo,1,oo_max); // operator series count: 4*w-1
       }

      Ref = I(0);
      for (int i = 0; i < w/2; ++i)
       {
           Ref = Ref(0) * Cell.swap_gate_no_sign(i, w-i-1);
           OpProgress(&oo,1,oo_max); // operator series count: 4*w+(w/2)-1
       }

      // if we could have tripartite symmetry, add operators for the sublattice magnetization
      //UnitCellMPO S_A, S_B, S_C;
      //not yet implemented ...

      UnitCellMPO Ry = I(0);
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
      OpProgress(&oo,1,oo_max); // operator series count: 4*w+(w/2)


      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO Hj, Hx_f, H_intra2;

      for (int i = 0; i < w; ++i)
      {
        // TIM - NN bonds
        Hj += Sz(0)[i] * Sz(0)[(i+1)%w]; // vertical bonds
        Hj += Sz(0)[i] * Sz(1)[i];       // -60 degree bonds
        Hj += Sz(0)[i] * Sz(1)[(i+1)%w]; // +60 degree bonds

        OpProgress(&oo,3,oo_max); // operator series count: 7*w+(w/2)

        // TIM - x-field terms
        Hx_f += Sx(0)[i];

        OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w/2)

        for (int j = i+1; j < w; ++j)
        {
          // Long-range Ising -- inter-cell interations only for alpha=2.0
          H_intra2 += std::pow( ( std::sin( PI/w ) / std::sin( (j-i)*PI/w ) ) , 2) * Sz(0)[i]*Sz(0)[j];
          OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)
        }
      }

      //Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J"] = sum_unit(Hj);
      Lattice["Hx_field"] = sum_unit(Hx_f);
      Lattice["H_LongRangeIsing_IntraCell2"] = sum_unit(H_intra2);

      OpProgress(&oo,3,oo_max); // operator series count: 8*w+(w*w/2)+3


      UnitCellMPO RyOld = I(0);
      for (int c = 0; c < w; ++c)
      {
         UnitCellMPO ThisR = I(0);
         // get the 'pivot' site/bond that we reflect about
         int const p1 = (w-c-1)/2;
         int const p2 = (w-c)/2;

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
         RyOld = RyOld * ThisR;
      }

      Lattice["RyOld"] = prod_unit_left_to_right(RyOld.MPO(), w*w);

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+4

      // Momentum operators in Y direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);
      Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w*w);

      OpProgress(&oo,2,oo_max); // operator series count: 8*w+(w*w/2)+6

      // for even size unit cell, add rotation by pi
      if (w%2 == 0)
      {
         UnitCellMPO TyPi = I(0);
         for (int i = 0; i < w/2; ++i)
         {
            TyPi = TyPi * Cell.swap_gate_no_sign(i, i+w/2);
         }
         Lattice["TyPi"] = prod_unit_left_to_right(TyPi.MPO(), w);
      }

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+7

      Lattice.func("TIM")(arg("J") = 0.0, arg("Gamma") = 0.0)
        = "J*H_J + Gamma*Hx_field";

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+8

      // a basic function for Haldane-Shastry model with Sz*Sz interations
      Lattice.func("HS")(arg("lambda") = 0.5, arg("p") = "0", arg("q") = "0")
                  = "exp(-lambda)*sum_string_inner( Sz(0)[p], exp(-lambda)*I(0), Sz(0)[q] )";

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+9

      Lattice.func("LongRangeIsing_IntraCell_YC4")(arg("alpha") = 2.0)
                  = "sum_unit( ( ( sin( pi/4 ) / sin( (1-0)*pi/4 ) )^alpha * Sz(0)[0]*Sz(0)[1] ) + "
                    "( ( sin( pi/4 ) / sin( (2-0)*pi/4 ) )^alpha * Sz(0)[0]*Sz(0)[2] ) + "
                    "( ( sin( pi/4 ) / sin( (3-0)*pi/4 ) )^alpha * Sz(0)[0]*Sz(0)[3] ) + "
                    "( ( sin( pi/4 ) / sin( (2-1)*pi/4 ) )^alpha * Sz(0)[1]*Sz(0)[2] ) + "
                    "( ( sin( pi/4 ) / sin( (3-1)*pi/4 ) )^alpha * Sz(0)[1]*Sz(0)[3] ) + "
                    "( ( sin( pi/4 ) / sin( (3-2)*pi/4 ) )^alpha * Sz(0)[2]*Sz(0)[3] ) )";

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+10

      Lattice.func("LongRangeIsing_IntraCell_YC6")(arg("alpha") = 2.0)
                  = "sum_unit( ( ( sin( pi/6 ) / sin( (1-0)*pi/6 ) )^alpha * Sz(0)[0]*Sz(0)[1] ) + "
                    "( ( sin( pi/6 ) / sin( (2-0)*pi/6 ) )^alpha * Sz(0)[0]*Sz(0)[2] ) + "
                    "( ( sin( pi/6 ) / sin( (3-0)*pi/6 ) )^alpha * Sz(0)[0]*Sz(0)[3] ) + "
                    "( ( sin( pi/6 ) / sin( (4-0)*pi/6 ) )^alpha * Sz(0)[0]*Sz(0)[4] ) + "
                    "( ( sin( pi/6 ) / sin( (5-0)*pi/6 ) )^alpha * Sz(0)[0]*Sz(0)[5] ) + "
                    "( ( sin( pi/6 ) / sin( (2-1)*pi/6 ) )^alpha * Sz(0)[1]*Sz(0)[2] ) + "
                    "( ( sin( pi/6 ) / sin( (3-1)*pi/6 ) )^alpha * Sz(0)[1]*Sz(0)[3] ) + "
                    "( ( sin( pi/6 ) / sin( (4-1)*pi/6 ) )^alpha * Sz(0)[1]*Sz(0)[4] ) + "
                    "( ( sin( pi/6 ) / sin( (5-1)*pi/6 ) )^alpha * Sz(0)[1]*Sz(0)[5] ) + "
                    "( ( sin( pi/6 ) / sin( (3-2)*pi/6 ) )^alpha * Sz(0)[2]*Sz(0)[3] ) + "
                    "( ( sin( pi/6 ) / sin( (4-2)*pi/6 ) )^alpha * Sz(0)[2]*Sz(0)[4] ) + "
                    "( ( sin( pi/6 ) / sin( (5-2)*pi/6 ) )^alpha * Sz(0)[2]*Sz(0)[5] ) + "
                    "( ( sin( pi/6 ) / sin( (4-3)*pi/6 ) )^alpha * Sz(0)[3]*Sz(0)[4] ) + "
                    "( ( sin( pi/6 ) / sin( (5-3)*pi/6 ) )^alpha * Sz(0)[3]*Sz(0)[5] ) + "
                    "( ( sin( pi/6 ) / sin( (5-4)*pi/6 ) )^alpha * Sz(0)[4]*Sz(0)[5] ) )";

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+11

      Lattice.func("LongRangeIsing_InterCell_YC4")(arg("alpha00")=0.0, arg("lambda00")=0.0, arg("alpha01")=0.0, arg("lambda01")=0.0, arg("alpha02")=0.0, arg("lambda02")=0.0, arg("alpha03")=0.0, arg("lambda03")=0.0)
                  = "alpha00*(HS{lambda=lambda00,p=0,q=0} + HS{lambda=lambda00,p=1,q=1} + HS{lambda=lambda00,p=2,q=2} + HS{lambda=lambda00,p=3,q=3}) + "
                    "alpha01*(HS{lambda=lambda01,p=0,q=1} + HS{lambda=lambda01,p=1,q=2} + HS{lambda=lambda01,p=2,q=3} + HS{lambda=lambda01,p=3,q=0}) + "
                    "alpha02*(HS{lambda=lambda02,p=0,q=2} + HS{lambda=lambda02,p=1,q=3} + HS{lambda=lambda02,p=2,q=0} + HS{lambda=lambda02,p=3,q=1}) + "
                    "alpha03*(HS{lambda=lambda03,p=0,q=3} + HS{lambda=lambda03,p=1,q=0} + HS{lambda=lambda03,p=2,q=1} + HS{lambda=lambda03,p=3,q=2})";

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+12

      Lattice.func("LongRangeIsing_InterCell_YC6_part1")(arg("alpha00")=0.0, arg("lambda00")=0.0, arg("alpha01")=0.0, arg("lambda01")=0.0, arg("alpha02")=0.0, arg("lambda02")=0.0, arg("alpha03")=0.0, arg("lambda03")=0.0, arg("alpha04")=0.0, arg("lambda04")=0.0, arg("alpha05")=0.0, arg("lambda05")=0.0)
                  = "alpha00*(HS{lambda=lambda00,p=0,q=0} + HS{lambda=lambda00,p=1,q=1} + HS{lambda=lambda00,p=2,q=2} + HS{lambda=lambda00,p=3,q=3}  + HS{lambda=lambda00,p=4,q=4} + HS{lambda=lambda00,p=5,q=5}) + "
                    "alpha01*(HS{lambda=lambda01,p=0,q=1} + HS{lambda=lambda01,p=1,q=2} + HS{lambda=lambda01,p=2,q=3} + HS{lambda=lambda01,p=3,q=4}  + HS{lambda=lambda01,p=4,q=5} + HS{lambda=lambda01,p=5,q=0}) + "
                    "alpha02*(HS{lambda=lambda02,p=0,q=2} + HS{lambda=lambda02,p=1,q=3} + HS{lambda=lambda02,p=2,q=4} + HS{lambda=lambda02,p=3,q=5}  + HS{lambda=lambda02,p=4,q=0} + HS{lambda=lambda02,p=5,q=1}) + "
                    "alpha03*(HS{lambda=lambda03,p=0,q=3} + HS{lambda=lambda03,p=1,q=4} + HS{lambda=lambda03,p=2,q=5} + HS{lambda=lambda03,p=3,q=0}  + HS{lambda=lambda03,p=4,q=1} + HS{lambda=lambda03,p=5,q=2}) + "
                    "alpha04*(HS{lambda=lambda04,p=0,q=4} + HS{lambda=lambda04,p=1,q=5} + HS{lambda=lambda04,p=2,q=0} + HS{lambda=lambda04,p=3,q=1}  + HS{lambda=lambda04,p=4,q=2} + HS{lambda=lambda04,p=5,q=3}) + "
                    "alpha05*(HS{lambda=lambda05,p=0,q=5} + HS{lambda=lambda05,p=1,q=0} + HS{lambda=lambda05,p=2,q=1} + HS{lambda=lambda05,p=3,q=2}  + HS{lambda=lambda05,p=4,q=3} + HS{lambda=lambda05,p=5,q=4})";

      OpProgress(&oo,1,oo_max); // operator series count: 8*w+(w*w/2)+13

      std::cout << "--> finished.\n";

      // save the lattice
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
