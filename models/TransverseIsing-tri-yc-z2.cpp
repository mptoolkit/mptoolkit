// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/TransverseIsing-tri-yc-z2.cpp
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

//
// YC configuration of a triangular lattice.
// The default unit-cell size is the width value.
//

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-z2.h"
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
	 std::cerr << "Constructs a triangular lattice in the YC configuration with a\n"
		   << "efficient way of numbering 1D chain. Operators are:\n"
		   << "H_J      - nearest neighbor spin exchange\n"
		   << "H_gamma  - transverse field\n"
                   << "S_z      - total spin on a leg of cylinder\n"
	    ;
         return 1;
      }

      LatticeSite Site = CreateZ2SpinSite(Spin);
      UnitCell Cell = repeat(Site, w);
      UnitCellOperator Sz(Cell, "Sz"), Sx(Cell, "Sx"), S_z(Cell, "S_z");

      // Add some operators on the unit-cell
      for (int i = 0; i < w; ++i)
      {
	 S_z += Sz[i];     // total spin on a leg of cylinder
      }

      // Now we construct the InfiniteLattice
      InfiniteLattice Lattice(Cell);        

      // Construct the Hamiltonian for a single unit-cell
      UnitCellMPO Hj, Hgamma;
      for (int i = 0; i < w; ++i)
      {
	 // Nearest neighbor bonds
	 // vertical bonds
	 Hj += Sz(0)[i] * Sz(0)[(i+1)%w];
	 // 60 degree bonds
	 Hj += Sz(0)[i] * Sz(1)[i];
	 Hj += Sz(0)[i] * Sz(1)[(i+1)%w];

	 // transverse field
	 Hgamma += Sx(0)[i];   
      }

      Lattice["H_J"] = sum_unit(Hj);
      Lattice["H_gamma"] = sum_unit(Hgamma);

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
