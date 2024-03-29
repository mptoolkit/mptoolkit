// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spin-tri-yc-su2-NonEff.cpp
//
// Copyright (C) 2014-2016 Seyed Saadatmand <s.saadatmand@uq.edu.au>
// Copyright (C) 2014-2024 Ian McCulloch <ian@qusim.net>
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

// YC configuration of a triangular lattice
//
// Example for a width-6 lattice (site numbers in brackets are periodic repeats in the vertical
// direction (ie top-left (0) is the same site as the bottom-left 0).
// Sites 6,7,8,9,10,11 are the second unit cell, eg 6 is (1)[0]
//
// Note that sites 5,6 are NOT nearest neighbor.
//
//     (18)
//  (6)  |
//(0) | 17
// | 11< |
// 5< | 16
// | 10< |
// 4< | 15
// |  9< |
// 3< | 14
// |  8< |
// 2< | 13
// |  7< |
// 1< | 12
// |  6< |
// 0< |(17)
// |(11)
//(5)
//

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;
using math_const::pi;

std::complex<double> phase(double theta)
{
   return std::exp(std::complex<double>(0.0, theta));
}


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
      int w = 4;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("width,w", prog_opt::value(&w), "width of the cylinder [default 4]")
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
         std::cerr << "Constructs a triangular lattice in the YC configuration with lattice vector (0,1)\n"
                   << "Operators:\n"
                   << "H_J1    - nearest neighbor spin exchange\n"
                   << "H_J2    - next-nearest neighbor spin exchange\n"
                   << "\nIf the lattice is tripartite (width is a multiple of 3) then we define sublattice\n"
                   << "spin operators on width*3 unit cells,\n"
                   << "S_A     - tripartite sublattice spin, including site S(0)[0]\n"
                   << "S_B     - tripartite sublattice spin, including site S(0)[1]\n"
                   << "S_C     - tripartite sublattice spin, including site S(0)[2]\n"

            ;
         return 1;
      }

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, w);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator S(Cell, "S"), StagS(Cell, "StagS");
      UnitCellOperator H1(Cell, "H1");
      UnitCellOperator H2(Cell, "H2");

      // Add some operators on the unit cell
      for (int i = 0; i < w; ++i)
      {
         S += S[i];     // total spin
      }

      // Staggered magnetization requires width to be even
      if (w%2 == 0)
      {
        for (int i = 0; i < w; ++i)
        {
           // staggered magnetization with alternating sublattices in Y-direction
           // (note: only one of the three possible formations).
           StagS += IntPow(-1,i) * S[i];
        }
       }
                  
      // if we have tripartite symmetry, add operators for the sublattice magnetization
      if (w%3 == 0)
      {
         UnitCellMPO S_A, S_B, S_C, S_120;
         for (int i = 0; i < w; i += 3)
         {
            S_A += S(0)[i]   + S(1)[(i+1)%w] + S(2)[(i+2)%w];
            S_B += S(0)[i+1] + S(1)[(i+2)%w] + S(2)[i%w];
            S_C += S(0)[i+2] + S(1)[i%w]     + S(2)[(i+1)%w];
            S_120 += S_A + phase(4*pi/3)*S_B + phase(-4*pi/3)*S_C;
         }
         Lattice["Sa"] = sum_unit(S_A, w*3);
         Lattice["Sb"] = sum_unit(S_B, w*3);
         Lattice["Sc"] = sum_unit(S_C, w*3);
         Lattice["S120"] = sum_unit(S_120, w*3);
      }

      // Construct the Hamiltonian for a single unit cell
      //      UnitCellMPO H1, H2;
      for (int i = 0; i < w; ++i)
      {
         // Nearest neighbor bonds
         H1 += inner(S(0)[i], S(0)[(i+1)%w]);      // vertical bonds
         H1 += inner(S(0)[i], S(1)[i]);            // 60 degree bonds - up
         H1 += inner(S(0)[i], S(1)[(i+w-1)%w]);    // 60 degree bonds - down

         // next-nearest neighbor bonds
         H2 += inner(S(0)[i], S(2)[(i+w-1)%w]);     // horizontal
         H2 += inner(S(0)[i], S(1)[(i+1)%w]);       // up-right
         H2 += inner(S(0)[i], S(1)[(i+w-2)%w]);     // down-right
      }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);

      // Lattice.func("H")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0)
      //   = "J1*H_J1 + J2*H_J2";  // a lattice function for an old set of calculations during 2014-15 

      Lattice.func("THM2")(arg("J2") = 0.0)
        = "H_J1 + J2*H_J2";

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
