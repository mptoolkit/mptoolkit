// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/spin-tri-2SiteUnitCell-su2.cpp
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2015 Seyed Saadatmand <s.saadatmand@uq.edu.au>
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
         std::cerr << "Constructs a triangular lattice in the C_wrapping=(-2,w) cylindrical configuration with a 2-site size unit-cell,\n"
                   << "and efficient way of numbering for 1D chain. Operators are:\n"
                   << "S_w     - total spin on a leg of cylinder\n"
                   << "H_J1    - nearest neighbor spin exchange\n"
                   << "H_J2    - next-nearest neighbor spin exchange\n"
            ;
         return 1;
      }

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site,2);
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator S(Cell, "S");

      // Add some operators on the unit-cell
      UnitCellMPO S_w;
      for (int i = 0; i < w; ++i)
      {
         S_w += S(i)[0];     // total spin on a leg of cylinder
      }

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2;

      for (int i = 0; i < 2; ++i)
       {
         // Nearest neighbor bonds
         // vertical bonds
         H1 += inner( S(0)[i] , S(1)[i] );
         // 60 degree bonds
         H1 += inner( S(0)[i] , S(i*w)[(i+1)%2] );      // up-right
         H1 += inner( S(0)[i] , S(i*w+1)[(i+1)%2] );    // down-right

         // next-nearest neighbor bonds
         H2 += inner( S(0)[i] , S(w+1)[i] );            // horizontal
         H2 += inner( S(0)[i] , S(i*w-1)[(i+1)%2] );    // up-right
         H2 += inner( S(0)[i] , S(i*w+2)[(i+1)%2] );    // down-right
       }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);

      Lattice.func("H")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0)
         = "J1*H_J1 + J2*H_J2";

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
