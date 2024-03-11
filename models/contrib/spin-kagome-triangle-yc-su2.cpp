// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spin-kagome-triangle-yc-su2.cpp
//
// Copyright (C) 2015-2018 Ian McCulloch <ian@qusim.net>
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

// YC configuration of a kagome-triangle lattice
//
// This adds an additional triangle to each triangle of the Kagome lattice.
// The width parameter is the number of Kagome hexagons around the width of the
// cylinder (this is different to the definition in spin-kagome-yc-su2.cpp).
// The unit cell size is 9*width.
//
// Example for width 2 lattice (site numbers in brackets are periodic repeats in
// the vertical direction; i.e. top-left (18) is the same site as the bottom-left 18.
// Sites 19-to-35 are the second unit cell, [1](0..17)
//
//
//  |
//(13)         18        (40)
//  |          |\        /|
//  |          | \      / |
//  |          | 19    24 |
//  |          | /|\  /|\ |
//  |          |/ | \/ | \|
//(-19)        20 | 23 | 26
//  |          |\ | /\ | /|
//  |          | \|/  \|/ |
//  |          | 21    25 |
//  |          | /      \ |
//  |          |/        \|
//  0          22        (45)
//  |\        /|          |\
//  | \      / |          |
//  |  1    6  |          |
//  | /|\  /|\ |          |
//  |/ | \/ | \|          |/
//  2  |  5 |  8         (47)
//  |\ | /\ | /|          |\
//  | \|/  \|/ |          |
//  |  3    7  |          |
//  | /      \ |          |
//  |/        \|          |/
//  4         27        (49)
//  |          |\        /|
//  |          | \      / |
//  |          | 28    33 |
//  |          | /|\  /|\ |
//  |          |/ | \/ | \|
//(-10)        29 | 32 | 35
//  |          |\ | /\ | /|
//  |          | \|/  \|/ |
//  |          |  30   34 |
//  |          | /      \ |
//  |          |/        \|
//  9          31        (36)
//  |\        /|          |\
//  | \      / |
//  | 10    15 |
//  | /|\  /|\ |
//  |/ | \/ | \|
//  11 | 14 | 17
//  |\ | /\ | /|
//  | \|/  \|/ |
//  | 12    16 |
//  | /      \ |
//  |/        \|
//  13        (18)
//  |
// \|
// (-1)
// /|
// (0)
//
// The 'a' sites are the sites within the additional triangles.  These are sites
// (1,2,3) and (6,7,8) of the unit cell.  These have an interaction J_aa.
// We also have a nearest-neighbour interaction J_ab that go between the
// 'a' site of the triangle, and the 'b' sites that make up the usual Kagome lattice.
// This is sites (0,1), (0,2), (1,5), (3,5), (2,4), (3,4), (5,6), (5,7),
// (6,*4), (8,*4), (7,**0), (8,**0),
// where '*' denotes the next unit to the right-above, and ** denotes
// the unit cell to the right-below.

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

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
         ("width,w", prog_opt::value(&w), "width of the cylinder, should be even [default 4]")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      // Descriptions of each operator
      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("SU(2) spin chain");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.author("Yan-Wei Dai", "daiyanwei1027@126.com");
      OpDescriptions.add_operators()
         ("H_aa"  , "nearest neighbor spin exchange between a-a sites")
         ("H_ab"  , "nearest neighbor spin exchange between a-b sites")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      unsigned u = 9*w;  // Hamiltonian unit-cell size

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, u);
      InfiniteLattice Lattice(&Cell);

      // Add some operators on the unit cell
      UnitCellOperator S(Cell, "S");
      for (int i = 0; i < u; ++i)
      {
         S += S[i];     // total spin on a leg of cylinder
      }

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO Haa, Hab;

      // a-a bonds within each triangle
      for (int i = 0; i < w; ++i)
      {
         // first triangle
         Haa += inner(S(0)[i*9+1], S(0)[i*9+2]);
         Haa += inner(S(0)[i*9+2], S(0)[i*9+3]);
         Haa += inner(S(0)[i*9+3], S(0)[i*9+1]);
         // second triangle
         Haa += inner(S(0)[i*9+6], S(0)[i*9+7]);
         Haa += inner(S(0)[i*9+7], S(0)[i*9+8]);
         Haa += inner(S(0)[i*9+8], S(0)[i*9+6]);
      }

      // a-b bonds between the triangles and kagome sites
      // Firstly the parts strictly within a unit cell
      for (int i = 0; i < w; ++i)
      {
         Hab += inner(S(0)[i*9+0], S(0)[i*9+1]);
         Hab += inner(S(0)[i*9+0], S(0)[i*9+2]);
         Hab += inner(S(0)[i*9+1], S(0)[i*9+5]);
         Hab += inner(S(0)[i*9+3], S(0)[i*9+5]);
         Hab += inner(S(0)[i*9+2], S(0)[i*9+4]);
         Hab += inner(S(0)[i*9+3], S(0)[i*9+4]);
         Hab += inner(S(0)[i*9+5], S(0)[i*9+6]);
         Hab += inner(S(0)[i*9+5], S(0)[i*9+7]);
      }
      // now the parts that cross unit cell boundaries
      for (int i = 0; i < w-1; ++i)
      {
         // get the starting point for the unit cells to the right-up and right-down
         int u = (i+w) % w;
         int l = (i+1) % w;
         u *= 9;
         l *= 9;

         Hab += inner(S(0)[i*9+6], S(1)[u+4]);
         Hab += inner(S(0)[i*9+8], S(1)[u+4]);

         Hab += inner(S(0)[i*9+7], S(1)[l]);
         Hab += inner(S(0)[i*9+8], S(1)[l]);
      }

      std::cout << "UnitCell size is:" << " " << u << std::endl;

      Lattice["H_aa"] = sum_unit(Haa);
      Lattice["H_ab"] = sum_unit(Hab);

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
