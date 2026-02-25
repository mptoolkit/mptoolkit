// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spin-kagome-yc.cpp
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2016 Seyed Saadatmand <s.saadatmand@uq.edu.au>
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

// Description: spin systems on kagome lattices with YC structure and efficient way of numbering;
// U(1)-symmetric.

// YC configuration of a kagome lattice.
// The default unit-cell size is '(3/2)*width' value, so the width should be even.
//
// Example for a 'width=6' lattice (site numbers in brackets are periodic repeats in
// the vertical direction; i.e. top-left (8) is the same site as the bottom-left 5).
// Sites 9-to-17 are the second unit cell, e.g. 9 is (1)[0].
//
//                (26)
//                 |
//        (17)     18
//         |       | >19
//(8)      9       20
// |       | >10 < |    ...
// 0       11      21
// | > 1 < |       | >22
// 2       12      23
// |       | >13 < |    ...
// 3       14      24
// | > 4 < |       | >25
// 5       15      26
// |       | >16 < |
// 6       17     (18)
// | > 7 < |
// 8      (9)
// |
//(0)


#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin.h"
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

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("Kagome lattice spin model");
      OpDescriptions.author("Jesse Osborne", "jesse.osborne@mpq.mpg.de");
      OpDescriptions.add_operators()
         ("H_J1"       , "nearest neighbor spin exchange")
         ("H_xx"       , "nearest neighbor xx spin coupling")
         ("H_yy"       , "nearest neighbor yy spin coupling")
         ("H_zz"       , "nearest neighbor zz spin coupling")
         ("H_x"        , "magnetic field in the x direction")
         ("H_y"        , "magnetic field in the y direction")
         ("H_z"        , "magnetic field in the z direction")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      CHECK(w%2 == 0)(w)("Width must be an even integer!");

      unsigned u = 3*(w/2);  // Hamiltonian unit-cell size

      LatticeSite Site = SpinSite(Spin);
      UnitCell Cell = repeat(Site, u);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz");

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H_xx, H_yy, H_zz,
                  H_xx2, H_yy2, H_zz2,
                  H_x, H_y, H_z;
      for (int i = 0; i < u; ++i)
      {
         // Nearest neighbor bonds

         // vertical bonds:
         if ( i%3 != 1 )
         {
            if ( i%3 == 0 )
            {
               H_xx += Sx(0)[i]*Sx(0)[i+2];
               H_yy += Sy(0)[i]*Sy(0)[i+2];
               H_zz += Sz(0)[i]*Sz(0)[i+2];
            }
            else
            {
               H_xx += Sx(0)[i]*Sx(0)[(i+1)%u];
               H_yy += Sy(0)[i]*Sy(0)[(i+1)%u];
               H_zz += Sz(0)[i]*Sz(0)[(i+1)%u];
            }
         }

         // 60 degree bonds:
         if ( i%3 == 0 )
         {
            H_xx += Sx(0)[i]*Sx(0)[i+1];
            H_yy += Sy(0)[i]*Sy(0)[i+1];
            H_zz += Sz(0)[i]*Sz(0)[i+1];
         }
         else if ( i%3 == 1 )
         {
            H_xx += Sx(0)[i]*Sx(1)[i+1];
            H_yy += Sy(0)[i]*Sy(1)[i+1];
            H_zz += Sz(0)[i]*Sz(1)[i+1];

            H_xx += Sx(0)[i]*Sx(1)[(i+2)%u];
            H_yy += Sy(0)[i]*Sy(1)[(i+2)%u];
            H_zz += Sz(0)[i]*Sz(1)[(i+2)%u];
         }
         else if ( i%3 == 2 )
         {
            H_xx += Sx(0)[i]*Sx(0)[i-1];
            H_yy += Sy(0)[i]*Sy(0)[i-1];
            H_zz += Sz(0)[i]*Sz(0)[i-1];
         }

         // Next-nearest neighbor bonds:
         // horizontal bonds:
         if ( i%3 != 1 )
         {
            if ( i%3 == 0 )
            {
               H_xx2 += Sx(0)[i]*Sx(1)[i+2];
               H_yy2 += Sy(0)[i]*Sy(1)[i+2];
               H_zz2 += Sz(0)[i]*Sz(1)[i+2];
            }
            else
            {
               H_xx2 += Sx(0)[i]*Sx(1)[(i+1)%u];
               H_yy2 += Sy(0)[i]*Sy(1)[(i+1)%u];
               H_zz2 += Sz(0)[i]*Sz(1)[(i+1)%u];
            }
         }

         // inclined bonds:
         if ( i%3 == 0 )
         {
            H_xx2 += Sx(0)[i]*Sx(0)[(i-2+u)%u];
            H_yy2 += Sy(0)[i]*Sy(0)[(i-2+u)%u];
            H_zz2 += Sz(0)[i]*Sz(0)[(i-2+u)%u];
         }
         else if ( i%3 == 1 )
         {
            H_xx2 += Sx(0)[i]*Sx(1)[i-1];
            H_yy2 += Sy(0)[i]*Sy(1)[i-1];
            H_zz2 += Sz(0)[i]*Sz(1)[i-1];

            H_xx2 += Sx(0)[i]*Sx(1)[(i+4)%u];
            H_yy2 += Sy(0)[i]*Sy(1)[(i+4)%u];
            H_zz2 += Sz(0)[i]*Sz(1)[(i+4)%u];
         }
         else if ( i%3 == 2 )
         {
            H_xx2 += Sx(0)[i]*Sx(0)[(i+2)%u];
            H_yy2 += Sy(0)[i]*Sy(0)[(i+2)%u];
            H_zz2 += Sz(0)[i]*Sz(0)[(i+2)%u];
         }

         H_x += Sx(0)[i];
         H_y += Sy(0)[i];
         H_z += Sz(0)[i];
      }

      Lattice["H_xx"] = sum_unit(H_xx);
      Lattice["H_yy"] = sum_unit(H_yy);
      Lattice["H_zz"] = sum_unit(H_zz);
      Lattice["H_J1"] = sum_unit(H_xx+H_yy+H_zz);
      Lattice["H_xx2"] = sum_unit(H_xx2);
      Lattice["H_yy2"] = sum_unit(H_yy2);
      Lattice["H_zz2"] = sum_unit(H_zz2);
      Lattice["H_J2"] = sum_unit(H_xx2+H_yy2+H_zz2);
      Lattice["H_x"] = sum_unit(H_x);
      Lattice["H_y"] = sum_unit(H_y);
      Lattice["H_z"] = sum_unit(H_z);

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
