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

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 3;
      double theta = 0.0;
      std::string FileName;

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
         ("THM_flux"                             , "J1-J2 Heisenebrg Hamiltonian on a triangular lattice with twisted BC in Y-direction as exp(i*theta)")
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

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell = repeat(Site, w);

      std::cout << "Building all Hamiltonian operators:\n";

      // Add some operators on the unit cell
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"), Trans(Cell, "Trans");

      for (int i = 0; i < w; ++i)
      {
	 Sp += Sp[i];     // total S+ on a leg of cylinder
	 Sm += Sm[i];     // total S- on a leg of cylinder
	 Sz += Sz[i];     // total Sz on a leg of cylinder
         std::cout << "... " << std::flush;
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
      {
         Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
      }

      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO Hz_ver, Hz1, H1, H1_flux, H2, H2_flux, H_intra;
      for (int i = 0; i < w; ++i)
      {
         // TIM - nearest neighbor bonds

         Hz_ver += Sz(0)[i]*Sz(0)[(i+1)%w];             // --> vertical bonds
         Hz1    += Sz(0)[i]*Sz(0)[(i+1)%w];
         Hz1    += Sz(0)[i]*Sz(1)[i];                   // --> +60 degree bonds
         Hz1    += Sz(0)[i]*Sz(1)[(i+1)%w];             // --> -60 degree bonds
         std::cout << "... " << std::flush;

	 // TIM - nearest neighbor bonds

	 // --> vertical bonds
	 H1 += Sz(0)[i]*Sz(0)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%w] + Sm(0)[i]*Sp(0)[(i+1)%w]);
         std::cout << ". " << std::flush;

         if ( (i+1)%w == 0 )
           H1_flux += Sz(0)[i]*Sz(0)[(i+1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(0)[(i+1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(0)[(i+1)%w]);
         else
           H1_flux += Sz(0)[i]*Sz(0)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%w] + Sm(0)[i]*Sp(0)[(i+1)%w]);
         std::cout << ". " << std::flush;

	 // --> 60 degree bonds
	 H1 += Sz(0)[i]*Sz(1)[i] + 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
	 H1 += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+1)%w] + Sm(0)[i]*Sp(1)[(i+1)%w]);
         std::cout << ".. " << std::flush; 

         H1_flux += Sz(0)[i]*Sz(1)[i] + 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
         if ( (i+1)%w == 0 )
           H1_flux += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(1)[(i+1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(1)[(i+1)%w]);
         else
           H1_flux += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+1)%w] + Sm(0)[i]*Sp(1)[(i+1)%w]);
         std::cout << ".. " << std::flush;       

	 // THM - next-nearest neighbor bonds
	 H2 += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
	 H2 += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+w-1)%w] + Sm(0)[i]*Sp(1)[(i+w-1)%w]);
	 H2 += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+2)%w] + Sm(0)[i]*Sp(1)[(i+2)%w]);
         std::cout << "... " << std::flush;

         if ( (i+1)%w == 0 )
           H2_flux += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(2)[(i+1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(2)[(i+1)%w]);
         else
           H2_flux += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
         if ( i == 0 )
           H2_flux += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(1)[(i+w-1)%w] + std::exp(-jj)*Sm(0)[i]*Sp(1)[(i+w-1)%w]);
         else
           H2_flux += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+w-1)%w] + Sm(0)[i]*Sp(1)[(i+w-1)%w]);
         if ( (i+1)%w == 0 )
           H2_flux += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (std::exp(jj)*Sp(0)[i]*Sm(1)[(i+2)%w] + std::exp(-jj)*Sm(0)[i]*Sp(1)[(i+2)%w]);
         else
           H2_flux += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+2)%w] + Sm(0)[i]*Sp(1)[(i+2)%w]);
         std::cout << "... " << std::flush;     

        for (int j = i+1; j < w; ++j)
        {
          // Long-range Ising - inter-cell interations
          H_intra += std::pow( ( std::sin( PI/w ) / std::sin( (j-i)*PI/w ) ) , 2) * Sz(0)[i]*Sz(0)[j];
          std::cout << ". " << std::flush;
        }
      }
       
      Lattice["H_ver"] = sum_unit(Hz_ver);
      Lattice["H_z1"] = sum_unit(Hz1);
      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);
      Lattice["H_J1_flux"] = sum_unit(H1_flux);
      Lattice["H_J2_flux"] = sum_unit(H2_flux);
      Lattice["H_LongRangeIsing_IntraCell"] = sum_unit(H_intra);
      std::cout << "... " << std::flush;

      // Momentum operators in Y direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);
      std::cout << ". " << std::flush;

      Lattice.func("THM")(arg("J1") = "cos(theta2)", arg("J2") = "sin(theta2)", arg("theta2") = "atan(alpha)", arg("alpha") = 0.0)
                  = "J1*H_J1 + J2*H_J2";
      std::cout << ". " << std::flush;

      Lattice.func("THM_flux")(arg("J1") = "cos(theta2)", arg("J2") = "sin(theta2)", arg("theta2") = "atan(alpha)", arg("alpha") = 0.0)
                  = "J1*H_J1_flux + J2*H_J2_flux";
      std::cout << ". " << std::flush;

      // a basic function for Haldane-Shastry model with Sz*Sz interations
      Lattice.func("HS")(arg("lambda") = 0.5, arg("i") = "0", arg("j") = "0")
                  = "exp(-lambda)*sum_string_inner( Sz(0)[i], exp(-lambda)*I(0), Sz(0)[j] )";
      std::cout << ". " << std::flush;

      /* Lattice.func("LongRangeIsing_NoInterCell_YC4")(arg("alpha00") = 1.0, arg("lambda00") = 0.5, arg("alpha01") = 1.0, arg("lambda01") = 0.5,  arg("alpha02") = 1.0, arg("lambda02") = 0.5)
                  = "alpha00*HS{lambda=lambda00,i=0,j=0} + alpha00*HS{lambda=lambda00,i=1,j=1} + alpha00*HS{lambda=lambda00,i=2,j=2} + alpha00*HS{lambda=lambda00,i=3,j=3} + alpha01*HS{lambda=lambda01,i=0,j=1} + alpha01*HS{lambda=lambda01,i=1,j=2} + alpha01*HS{lambda=lambda01,i=2,j=3} + alpha01*HS{lambda=lambda01,i=3,j=1} + alpha02*HS{lambda=lambda02,i=0,j=2} + alpha02*HS{lambda=lambda02,i=1,j=3}"; */

      /* Lattice.func("LongRangeIsing_InterCell_YC4")(arg("alpha00")=0.0, arg("lambda00")=0.0, arg("alpha01")=0.0, arg("lambda01")=0.0, arg("alpha02")=0.0, arg("lambda02")=0.0, arg("alpha03")=0.0, arg("lambda03")=0.0)
                  = "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[0], exp(-lambda00)*I(0), Sz(0)[0] ) ) + alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[1], exp(-lambda00)*I(0), Sz(0)[1] ) ) + alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[2], exp(-lambda00)*I(0), Sz(0)[2] ) ) + alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[3], exp(-lambda00)*I(0), Sz(0)[3] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[0], exp(-lambda01)*I(0), Sz(0)[1] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[1], exp(-lambda01)*I(0), Sz(0)[2] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[2], exp(-lambda01)*I(0), Sz(0)[3] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[3], exp(-lambda01)*I(0), Sz(0)[0] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[0], exp(-lambda02)*I(0), Sz(0)[2] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[1], exp(-lambda02)*I(0), Sz(0)[3] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[2], exp(-lambda02)*I(0), Sz(0)[0] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[3], exp(-lambda02)*I(0), Sz(0)[1] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[0], exp(-lambda03)*I(0), Sz(0)[3] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[1], exp(-lambda03)*I(0), Sz(0)[0] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[2], exp(-lambda03)*I(0), Sz(0)[1] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[3], exp(-lambda03)*I(0), Sz(0)[2] ) )"; */

      Lattice.func("LongRangeIsing_InterCell_YC4_test00")(arg("alpha00")=0.0, arg("lambda00")=0.0)
                  = "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[0], exp(-lambda00)*I(0), Sz(0)[0] ) )";
      std::cout << ". " << std::flush;

      Lattice.func("LongRangeIsing_InterCell_YC4_test01")(arg("alpha01")=0.0, arg("lambda01")=0.0)
                  = "alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[0], exp(-lambda01)*I(0), Sz(0)[1] ) )";
      std::cout << ". " << std::flush;

      Lattice.func("LongRangeIsing_InterCell_YC4_part1")(arg("alpha00")=0.0, arg("lambda00")=0.0, arg("alpha01")=0.0, arg("lambda01")=0.0)
                  = "alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[0], exp(-lambda00)*I(0), Sz(0)[0] ) ) + alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[1], exp(-lambda00)*I(0), Sz(0)[1] ) ) + alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[2], exp(-lambda00)*I(0), Sz(0)[2] ) ) + alpha00*( exp(-lambda00)*sum_string_inner( Sz(0)[3], exp(-lambda00)*I(0), Sz(0)[3] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[0], exp(-lambda01)*I(0), Sz(0)[1] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[1], exp(-lambda01)*I(0), Sz(0)[2] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[2], exp(-lambda01)*I(0), Sz(0)[3] ) ) + alpha01*( exp(-lambda01)*sum_string_inner( Sz(0)[3], exp(-lambda01)*I(0), Sz(0)[0] ) )";
      std::cout << ". " << std::flush;

      Lattice.func("LongRangeIsing_InterCell_YC4_part2")(arg("alpha02")=0.0, arg("lambda02")=0.0, arg("alpha03")=0.0, arg("lambda03")=0.0)
                  = "alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[0], exp(-lambda02)*I(0), Sz(0)[2] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[1], exp(-lambda02)*I(0), Sz(0)[3] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[2], exp(-lambda02)*I(0), Sz(0)[0] ) ) + alpha02*( exp(-lambda02)*sum_string_inner( Sz(0)[3], exp(-lambda02)*I(0), Sz(0)[1] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[0], exp(-lambda03)*I(0), Sz(0)[3] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[1], exp(-lambda03)*I(0), Sz(0)[0] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[2], exp(-lambda03)*I(0), Sz(0)[1] ) ) + alpha03*( exp(-lambda03)*sum_string_inner( Sz(0)[3], exp(-lambda03)*I(0), Sz(0)[2] ) )";
      std::cout << ". " << std::flush;

      std::cout << ">>> finished.\n";

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
