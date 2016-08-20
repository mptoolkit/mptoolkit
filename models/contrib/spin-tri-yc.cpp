// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-tri-yc.cpp
//
// Copyright (C) 2015,2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2015,2016 Seyed N. Saadatmand <s.saadatmand@uq.edu.au>
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
#include "models/spin.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>


namespace prog_opt = boost::program_options;


int IntPow(int x, int p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * IntPow(x, p-1);
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
                   << "not yet implemented ...\n\n"
            ;
         return 1;
      }

      const double PI = 4.0*std::atan(1.0);
      //const std::complex<double> jj(0.0,theta*PI);
      int oo = 0;
      int oo_max = 11*w+(w*w/2)+14;

      LatticeSite Site = SpinSite(Spin);
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

         oo+=3;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 3*w
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
       {
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);

           oo++;
           std::printf("\33[2K\r");
           std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 4*w-1
       }

      Ref = I(0);
      for (int i = 0; i < w/2; ++i)
       {
           Ref = Ref(0) * Cell.swap_gate_no_sign(i, w-i-1);

           oo++;
           std::printf("\33[2K\r");
           std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 4*w+(w/2)-1
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

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 4*w+(w/2)

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO Hx_f, H_intra, H1, H2;

      for (int i = 0; i < w; ++i)
      {

        // TIM - x-field terms
        Hx_f += Sx(0)[i];

        oo++;
        std::printf("\33[2K\r");
        std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 5*w+(w/2)

        // THM - nearest neighbor bonds

        // --> vertical bonds
        H1 += Sx(0)[i]*Sx(0)[(i+1)%w] + Sy(0)[i]*Sy(0)[(i+1)%w] + Sz(0)[i]*Sz(0)[(i+1)%w];

        // --> 60 degree bonds
        H1 += Sx(0)[i]*Sx(1)[i] + Sy(0)[i]*Sy(1)[i] + Sz(0)[i]*Sz(1)[i];
        H1 += Sx(0)[i]*Sx(1)[(i+1)%w] + Sy(0)[i]*Sy(1)[(i+1)%w] + Sz(0)[i]*Sz(1)[(i+1)%w];

        oo++;
        std::printf("\33[2K\r");
        std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 8*w+(w/2)

        // THM - next-nearest neighbor bonds
        H2 += Sx(0)[i]*Sx(2)[(i+1)%w] + Sy(0)[i]*Sy(2)[(i+1)%w] + Sz(0)[i]*Sz(2)[(i+1)%w];
        H2 += Sx(0)[i]*Sx(1)[(i+w-1)%w] + Sy(0)[i]*Sy(1)[(i+w-1)%w] + Sz(0)[i]*Sz(1)[(i+w-1)%w];
        H2 += Sx(0)[i]*Sx(1)[(i+2)%w] + Sy(0)[i]*Sy(1)[(i+2)%w] + Sz(0)[i]*Sz(1)[(i+2)%w];

        oo++;
        std::printf("\33[2K\r");
        std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w/2)

        for (int j = i+1; j < w; ++j)
        {
          // Long-range Ising - inter-cell interations
          H_intra += std::pow( ( std::sin( PI/w ) / std::sin( (j-i)*PI/w ) ) , 2) * Sz(0)[i]*Sz(0)[j];

          oo++;
          std::printf("\33[2K\r");
          std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)
        }
      }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);
      Lattice["Hx_field"] = sum_unit(Hx_f);
      Lattice["H_LongRangeIsing_IntraCell"] = sum_unit(H_intra);

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+4


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

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+5

      // Momentum operators in Y direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);
      Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w*w);

      oo+=2;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+7


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

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+8

      Lattice.func("THM2")(arg("J2") = 0.0)
              = "H_J1 + J2*H_J2";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+9

      // a basic function for Haldane-Shastry model with Sz*Sz interations
      Lattice.func("HS")(arg("lambda") = 0.5, arg("i") = "0", arg("j") = "0")
                  = "exp(-lambda)*sum_string_inner( Sz(0)[i], exp(-lambda)*I(0), Sz(0)[j] )";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+10

      Lattice.func("LongRangeIsing_InterCell_YC4_part1")(arg("alpha00")=0.0, arg("lambda00")=0.0, arg("alpha01")=0.0, arg("lambda01")=0.0)
                  = "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[0], exp(-lambda00)*I(0), Sz(0)[0] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[1], exp(-lambda00)*I(0), Sz(0)[1] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[2], exp(-lambda00)*I(0), Sz(0)[2] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[3], exp(-lambda00)*I(0), Sz(0)[3] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[0], exp(-lambda01)*I(0), Sz(0)[1] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[1], exp(-lambda01)*I(0), Sz(0)[2] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[2], exp(-lambda01)*I(0), Sz(0)[3] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[3], exp(-lambda01)*I(0), Sz(0)[0] ) )   ";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+11

      Lattice.func("LongRangeIsing_InterCell_YC4_part2")(arg("alpha02")=0.0, arg("lambda02")=0.0, arg("alpha03")=0.0, arg("lambda03")=0.0)
                  = "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[0], exp(-lambda02)*I(0), Sz(0)[2] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[1], exp(-lambda02)*I(0), Sz(0)[3] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[2], exp(-lambda02)*I(0), Sz(0)[0] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[3], exp(-lambda02)*I(0), Sz(0)[1] ) ) + "
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[0], exp(-lambda03)*I(0), Sz(0)[3] ) ) + "
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[1], exp(-lambda03)*I(0), Sz(0)[0] ) ) + "
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[2], exp(-lambda03)*I(0), Sz(0)[1] ) ) + "
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[3], exp(-lambda03)*I(0), Sz(0)[2] ) )   ";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+12

       Lattice.func("LongRangeIsing_InterCell_YC6_part1")(arg("alpha00")=0.0, arg("lambda00")=0.0, arg("alpha01")=0.0, arg("lambda01")=0.0, arg("alpha02")=0.0, arg("lambda02")=0.0)
                  = "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[0], exp(-lambda00)*I(0), Sz(0)[0] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[1], exp(-lambda00)*I(0), Sz(0)[1] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[2], exp(-lambda00)*I(0), Sz(0)[2] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[3], exp(-lambda00)*I(0), Sz(0)[3] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[4], exp(-lambda00)*I(0), Sz(0)[4] ) ) + "
                    "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[5], exp(-lambda00)*I(0), Sz(0)[5] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[0], exp(-lambda01)*I(0), Sz(0)[1] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[1], exp(-lambda01)*I(0), Sz(0)[2] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[2], exp(-lambda01)*I(0), Sz(0)[3] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[3], exp(-lambda01)*I(0), Sz(0)[4] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[4], exp(-lambda01)*I(0), Sz(0)[5] ) ) + "
                    "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[5], exp(-lambda01)*I(0), Sz(0)[0] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[0], exp(-lambda02)*I(0), Sz(0)[2] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[1], exp(-lambda02)*I(0), Sz(0)[3] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[2], exp(-lambda02)*I(0), Sz(0)[4] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[3], exp(-lambda02)*I(0), Sz(0)[5] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[4], exp(-lambda02)*I(0), Sz(0)[0] ) ) + "
                    "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[5], exp(-lambda02)*I(0), Sz(0)[1] ) )   ";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+13

      Lattice.func("LongRangeIsing_InterCell_YC6_part2")(arg("alpha03")=0.0, arg("lambda03")=0.0, arg("alpha04")=0.0, arg("lambda04")=0.0, arg("alpha05")=0.0, arg("lambda05")=0.0)
                  = "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[0], exp(-lambda03)*I(0), Sz(0)[3] ) ) +"
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[1], exp(-lambda03)*I(0), Sz(0)[4] ) ) +"
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[2], exp(-lambda03)*I(0), Sz(0)[5] ) ) +"
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[3], exp(-lambda03)*I(0), Sz(0)[0] ) ) +"
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[4], exp(-lambda03)*I(0), Sz(0)[1] ) ) +"
                    "alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[5], exp(-lambda03)*I(0), Sz(0)[2] ) ) +"
                    "alpha04*( exp(-lambda04)*sum_string_inner( Sz(0)[0], exp(-lambda04)*I(0), Sz(0)[4] ) ) +"
                    "alpha04*( exp(-lambda04)*sum_string_inner( Sz(0)[1], exp(-lambda04)*I(0), Sz(0)[5] ) ) +"
                    "alpha04*( exp(-lambda04)*sum_string_inner( Sz(0)[2], exp(-lambda04)*I(0), Sz(0)[0] ) ) +"
                    "alpha04*( exp(-lambda04)*sum_string_inner( Sz(0)[3], exp(-lambda04)*I(0), Sz(0)[1] ) ) +"
                    "alpha04*( exp(-lambda04)*sum_string_inner( Sz(0)[4], exp(-lambda04)*I(0), Sz(0)[2] ) ) +"
                    "alpha04*( exp(-lambda04)*sum_string_inner( Sz(0)[5], exp(-lambda04)*I(0), Sz(0)[3] ) ) +"
                    "alpha05*( exp(-lambda05)*sum_string_inner( Sz(0)[0], exp(-lambda05)*I(0), Sz(0)[5] ) ) +"
                    "alpha05*( exp(-lambda05)*sum_string_inner( Sz(0)[1], exp(-lambda05)*I(0), Sz(0)[0] ) ) +"
                    "alpha05*( exp(-lambda05)*sum_string_inner( Sz(0)[2], exp(-lambda05)*I(0), Sz(0)[1] ) ) +"
                    "alpha05*( exp(-lambda05)*sum_string_inner( Sz(0)[3], exp(-lambda05)*I(0), Sz(0)[2] ) ) +"
                    "alpha05*( exp(-lambda05)*sum_string_inner( Sz(0)[4], exp(-lambda05)*I(0), Sz(0)[3] ) ) +"
                    "alpha05*( exp(-lambda05)*sum_string_inner( Sz(0)[5], exp(-lambda05)*I(0), Sz(0)[4] ) )  ";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w+(w*w/2)+14

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
