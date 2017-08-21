// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-tri-yc-u1.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
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
// |  8< |
// 2< | 15
// |  9< |
// 3< | 16
// | 10< |
// 4< | 17
// | 11< |
// 5< | (12)
// | (6)
//(0)
//

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-u1.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

/* void ClearScreen()
{
    std::cout << std::string( 100, '\n' );
} */

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 3;
      double theta = 0.0;
      std::string FileName;
      bool NoReflect = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("width,w", prog_opt::value(&w), "width of the cylinder [default 3]")
         ("theta,t", prog_opt::value(&theta), "flux phase to twist boundary condition in Y-direction (in unit of PI) [default 0.0]")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.add_operators()
         ("H_ver"                  , "nearest neighbor Ising spin exchange in Y-direction")
         ("H_z1"                   , "nearest neighbor Ising spin exchange")
         ("H_J1"                   , "nearest neighbor Heisenberg spin exchange")
         ("H_J2"                   , "next-nearest neighbor Heisenberg spin exchange")
         ("H_J1_flux"              , "nearest neighbor Heisenberg spin exchange with flux")
         ("H_J2_flux"              , "next-nearest neighbor Heisenberg spin exchange with flux")
         ("H_LongRangeIsing_intra" , "intra-cell interactions of long-range Ising model")
         ;

      OpDescriptions.add_functions()
         ("THM_flux"                             , "J1-J2 Heiseneberg Hamiltonian on a triangular lattice with twisted BC in Y-direction as exp(i*theta)")
         ("HS"                                   , "Haldane-Shastry Hamiltonian with Sz*Sz interactions, parametized by 'lambda' (exponential decay as exp(-lambda*r))")
         ("LongRangeIsing_InterCell_YC4_part1"   , "long-range Ising model Hamiltonian on a 4-leg YC structure, parametized by 'alpha0j's and 'lambda0j's | PART 1")
         ("LongRangeIsing_InterCell_YC4_part2"   , "long-range Ising model Hamiltonian on a 4-leg YC structure, parametized by 'alpha0j's and 'lambda0j's | PART 2")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Constructs a triangular lattice in the YC configuration with lattice vector (0,1).\n";
         std::cerr << OpDescriptions << '\n';
         std::cerr << "only if the lattice could be potentially tripartite (width is a multiple of 3):\n"
                   << "define sublattice spin operators on 'width*3' unit-cells.\n";
            ;
         return 1;
      }

      const double PI = 4.0*std::atan(1.0);
      const std::complex<double> jj(0.0,theta*PI);
      int oo = 0;
      int oo_max = 21*w+(w*(w-1)/2)+18;

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell = repeat(Site, w);
      InfiniteLattice Lattice(&Cell);

      std::cout << "Building all Hamiltonian operators:\n";

      // Add some operators on the unit cell
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"), Trans(Cell, "Trans");

      for (int i = 0; i < w; ++i)
      {
         Sp += Sp[i];     // total S+ on a leg of cylinder
         Sm += Sm[i];     // total S- on a leg of cylinder
         Sz += Sz[i];     // total Sz on a leg of cylinder

         oo+=2;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 3*w
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
      {
         Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);

         oo++;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 4*w
      }

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
      }


      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO Hz_v, Hz1, H1, H1_flux, H2, H2_flux, H_intra2;

      for (int i = 0; i < w; ++i)
      {
         // TIM - x-field terms
         // Hx_f += 0.5*(Sp(0)[i]+Sm(0)[i]);

         oo++;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 5*w


         // TIM - nearest neighbor bonds

         Hz_v += Sz(0)[i]*Sz(0)[(i+1)%w];             // --> vertical bonds for 'Hz_ver'
         Hz1  += Sz(0)[i]*Sz(0)[(i+1)%w];             // --> vertical bonds for 'Hz1'
         Hz1  += Sz(0)[i]*Sz(1)[i];                   // --> +60 degree bonds
         Hz1  += Sz(0)[i]*Sz(1)[(i+1)%w];             // --> -60 degree bonds

         oo+=4;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 9*w


         // THM - nearest neighbor bonds

         // --> vertical bonds
         H1 += Sz(0)[i]*Sz(0)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%w] + Sm(0)[i]*Sp(0)[(i+1)%w]);

         oo++;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 10*w

         if ( (i+1)%w == 0 )
           H1_flux += Sz(0)[i]*Sz(0)[(i+1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(0)[(i+1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(0)[(i+1)%w]);
         else
           H1_flux += Sz(0)[i]*Sz(0)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%w] + Sm(0)[i]*Sp(0)[(i+1)%w]);

         oo++;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 11*w

         // --> 60 degree bonds
         H1 += Sz(0)[i]*Sz(1)[i] + 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
         H1 += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+1)%w] + Sm(0)[i]*Sp(1)[(i+1)%w]);

         oo+=2;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 13*w

         H1_flux += Sz(0)[i]*Sz(1)[i] + 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
         if ( (i+1)%w == 0 )
           H1_flux += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(1)[(i+1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(1)[(i+1)%w]);
         else
           H1_flux += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+1)%w] + Sm(0)[i]*Sp(1)[(i+1)%w]);

         oo+=2;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 15*w

         // THM - next-nearest neighbor bonds
         H2 += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
         H2 += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+w-1)%w] + Sm(0)[i]*Sp(1)[(i+w-1)%w]);
         H2 += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+2)%w] + Sm(0)[i]*Sp(1)[(i+2)%w]);

         oo+=3;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 18*w

         if ( (i+1)%w == 0 )
         {
           H2_flux += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(2)[(i+1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(2)[(i+1)%w]);
           H2_flux += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(1)[(i+2)%w] + std::exp(-jj)*Sm(0)[i]*Sp(1)[(i+2)%w]);
         }
         else
         {
           H2_flux += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
           H2_flux += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+2)%w] + Sm(0)[i]*Sp(1)[(i+2)%w]);
         }

         oo+=2;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 20*w

         if ( i == 0 )
           H2_flux += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(1)[(i+w-1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(1)[(i+w-1)%w]);
         else
           H2_flux += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+w-1)%w] + Sm(0)[i]*Sp(1)[(i+w-1)%w]);

         oo++;
         std::printf("\33[2K\r");
         std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w

        for (int j = i+1; j < w; ++j)
        {
          // Long-range Ising - inter-cell interations
          H_intra2 += std::pow( ( std::sin( PI/w ) / std::sin( (j-i)*PI/w ) ) , 2) * Sz(0)[i]*Sz(0)[j];

          oo++;
          std::printf("\33[2K\r");
          std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)
        }
      }

      Lattice["Hz_ver"] = sum_unit(Hz_v);
      Lattice["H_z1"] = sum_unit(Hz1);
      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);
      Lattice["H_J1_flux"] = sum_unit(H1_flux);
      Lattice["H_J2_flux"] = sum_unit(H2_flux);
      Lattice["H_LongRangeIsing_IntraCell2"] = sum_unit(H_intra2);

      oo+=8;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+8

      // Momentum operators in Y direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);

      // Reflection about Y
      Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w*w);

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+9

      Lattice.func("THM")(arg("J1") = "cos(theta2)", arg("J2") = "sin(theta2)", arg("theta2") = "atan(alpha)", arg("alpha") = 0.0)
                  = "J1*H_J1 + J2*H_J2"; // an old lattice function, used in few projects in 2014-15.

      Lattice.func("THM2")(arg("J2") = 0.0)
              = "H_J1 + J2*H_J2";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+11

      Lattice.func("THM_flux")(arg("J1") = "cos(theta2)", arg("J2") = "sin(theta2)", arg("theta2") = "atan(alpha)", arg("alpha") = 0.0)
                  = "J1*H_J1_flux + J2*H_J2_flux";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+12

      // a basic function for Haldane-Shastry model with Sz*Sz interations
      Lattice.func("HS")(arg("lambda") = 0.5, arg("j") = "0", arg("k") = "0")
                  = "exp(-lambda)*sum_string_inner( Sz(0)[j], exp(-lambda)*I(0), Sz(0)[k] )";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+13

      Lattice.func("TestFunc00")(arg("alpha00")=0.0, arg("lambda00")=0.0)
                  = "sin(2.5)*alpha00*HS{lambda=-lambda00,j=1,k=2}";

      Lattice.func("LongRangeIsing_IntraCell_YC4")(arg("alpha") = 2.0)
                  = "sum_unit( ( ( sin( pi/4 ) / sin( (1-0)*pi/4 ) )^alpha * Sz(0)[0]*Sz(0)[1] ) + "
                    "( ( sin( pi/4 ) / sin( (2-0)*pi/4 ) )^alpha * Sz(0)[0]*Sz(0)[2] ) + "
                    "( ( sin( pi/4 ) / sin( (3-0)*pi/4 ) )^alpha * Sz(0)[0]*Sz(0)[3] ) + "
                    "( ( sin( pi/4 ) / sin( (2-1)*pi/4 ) )^alpha * Sz(0)[1]*Sz(0)[2] ) + "
                    "( ( sin( pi/4 ) / sin( (3-1)*pi/4 ) )^alpha * Sz(0)[1]*Sz(0)[3] ) + "
                    "( ( sin( pi/4 ) / sin( (3-2)*pi/4 ) )^alpha * Sz(0)[2]*Sz(0)[3] ) )";

      oo++;
      std::printf("\33[2K\r");
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+14

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
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+15

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
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+16

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
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+17

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
      std::cout << "working... %" << (100*oo)/oo_max << std::flush; // operator series count: 21*w+(w*(w-1)/2)+18

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
