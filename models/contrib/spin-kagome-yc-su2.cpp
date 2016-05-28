// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-kagome-yc-su2.cpp
//
// Copyright (C) 2015,2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2015,2016 Seyed N. Saadatmand <s.saadatmand@uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
// Descriptin: spin systems on kagome lattices with YC structure and efficient way of numbering; SU(2)-symmetric. <spin-kagome-yc-su2.cpp>
// Authors: Ian P. McCulloch and Seyed N. Saadatmand
// Contact: s.saadatmand@uq.edu.au

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
      
      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Constructs a kagome lattice in the YC configuration with an\n"
		   << "efficient way of numbering of 1D chain. The default unit-cell size\n" 
                   << "is '(3/2)*width' value. SU(2)-symmetry is implied.\n"
                   << "Operators:\n"
		   << "H_J1    - nearest neighbor spin exchange\n"
		   << "H_J2    - next-nearest neighbor spin exchange\n"
                   << "Functions:\n"
                   << "H( J1 = NN coupling strength, J2 = NNN coupling strength, theta = atan(J2/J1)\n"
                   << "  \"radians\", alpha = J2/J1 )\n"   
	    ;
         return 1;
      }

      CHECK(w%2 == 0)(w)("Width must be an even integer!");

      unsigned u = 3*(w/2);  // Hamiltonian unit-cell size

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, u);

      // Add some operators on the unit cell
      UnitCellOperator S(Cell, "S");
      for (int i = 0; i < u; ++i)
      {
	 S += S[i];     // total spin on a leg of cylinder
      }

      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      unsigned num_bonds_j1 = 0;
      unsigned num_bonds_j2 = 0;

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2;
      for (int i = 0; i < u; ++i)
      {
	 // Nearest neighbor bonds

	 // vertical bonds:
	 if ( i%3 != 1 ) 
         {
           if ( i%3 == 0 )
             H1 += inner(S(0)[i], S(0)[i+2]);
           else
             H1 += inner(S(0)[i], S(0)[(i+1)%u]);
           ++num_bonds_j1;
         }
  
	 // 60 degree bonds:
	 if ( i%3 == 0 )
           {
             H1 += inner(S(0)[i], S(0)[i+1]);
             ++num_bonds_j1;
           }
         else if ( i%3 == 1 )
           {
             H1 += inner(S(0)[i], S(1)[i+1]);   
	     H1 += inner(S(0)[i], S(1)[(i+2)%u]);
             num_bonds_j1 += 2;
           }
         else if ( i%3 == 2 )
           {
             H1 += inner(S(0)[i], S(0)[i-1]);
             ++num_bonds_j1;
           }
 
	 // Next-nearest neighbor bonds
	 
	 // horizental bonds:
	 if ( i%3 != 1 )
         {
           if ( i%3 == 0 )
             H2 += inner(S(0)[i], S(1)[i+2]);
           else   
             H2 += inner(S(0)[i], S(1)[(i+1)%u]);
           ++num_bonds_j2;
         }

         // inclined bonds:
         if ( i%3 == 0 )
           {
             H2 += inner(S(0)[i], S(0)[(i-2+u)%u]);
             ++num_bonds_j2;
           }
         else if ( i%3 == 1 )
           {
             H2 += inner(S(0)[i], S(1)[i-1]);
             H2 += inner(S(0)[i], S(1)[(i+4)%u]);
             num_bonds_j2 += 2;
           }
         else if ( i%3 == 2 )
           {
             H2 += inner(S(0)[i], S(0)[(i+2)%u]);
             ++num_bonds_j2;
           }
      }

      std::cout << "UnitCell size is:" << " " << u << std::endl;  
      std::cout << "The number of J1 bonds per unit-cell:" << " " << num_bonds_j1 << ", and the number of J2 bonds per unit-cell: " << num_bonds_j2 << std::endl;

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
