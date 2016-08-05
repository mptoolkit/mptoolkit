// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-tri-xc-su2.cpp
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

// Description: spin systems on triangular lattices XC SU(2)-symmetric.
// XC configuration of a triangular lattice.
// The width W measures the number of total rows, which is twice the number of sites
// in a Y column.
// The unit-cell is two columns, of total size W
//
// Example for W=8:
//
//  (4)-(12)
//  / \ / 
//(0)-(8)--
//  \ / \ /
//   7---15
//  / \ / 
// 3---11--
//  \ / \ /
//   6---14
//  / \ / 
// 2---10--
//  \ / \ /
//   5---13
//  / \ / 
// 1---9---
//  \ / \ /
//   4---12
//  / \ / 
// 0---8---
//  \ / \ /
//  (7)-(15)
//  / \ / 
//(3)-(11)--

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
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
      int w = 6;
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
      
      OperatorDescriptions OpDescriptions;
      OpDescriptions.add_operators()
	 ("H_J1",     "nearest neighbor spin exchange")
	 ("H_J2",     "next-nearest neighbor spin exchange")
	 ("H_Jcell",  "zig-zag cylinder coupling")
	 ("Ty"  ,     "Translation in Y direction")
	 ("TyPi",     "Translation by pi in Y-direction (only if w is divisible by 4)")
	 ("Ry"  ,     "Reflection about the X-axis")
	 ;

      OpDescriptions.add_functions()
         ("THM",        "J1-J2 Heisenebrg Hamiltonian on a triangular lattice | inputs: J1 = cos(theta), J2 = sin(theta), theta = atan(alpha), alpha = 0.0")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << OpDescriptions << '\n';
         return 1;
      }

      if (w%2 != 0)
      {
	 std::cerr << "fatal: width must be even.\n";
	 return 1;
      }

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, w);
      UnitCellOperator S(Cell, "S"), StagS(Cell, "StagS");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      // w/2, for convenience
      int const w2 = w/2;

      // Add some operators on the unit-cell
   
      for (int i = 0; i < w; ++i)
      {
	 S += S[i];                       // total spin on a unit cell
         StagS += IntPow(-1,i/w2) * S[i]; // staggered magnetization in Y-direction (note: only one of the three possible formations).
      }

      Trans = I(0);
      for (int i = 0; i < w2-1; ++i)
       {
           //T *= 0.5*( 0.25*inner(S[i],S[i+1]) + 1 );
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
           Trans = Trans(0) * Cell.swap_gate_no_sign(i+w2, i+w2+1);
       }

       // to test existence of bipartite/tripartite symmetry, add operators for the staggered/sublattice magnetization:
       // NOTE: XC structure does NOT need a 3*W unit-cell for sublattice magnetization operators. 
       UnitCellMPO S_A, S_B, S_C;
       UnitCellMPO S_stag_n60;             

       for (int i = 0; i < w; i += 3)
         {
            S_A += S(0)[i];
            if ( (i+1)<w ) S_B += S(0)[i+1];
            if ( (i+2)<w ) S_C += S(0)[i+2];
            S_stag_n60 += IntPow(-1,i+(i/w2))*S(0)[i] + IntPow(-1,i+(i/w2)+1)*S(1)[i];    // staggered magnetization order parameter with FM stripes in -60^degree direction.  
         }
              
      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2, H_Jcell;
    
      for (int i = 0; i < w2; ++i)
      {
	 // nearest neighbor bonds, first column
	 H1 += inner(S(0)[i], S(1)[i]);               // horizontal
	 H1 += inner(S(0)[i], S(0)[i+w2]);            // up-right
	 H1 += inner(S(0)[i], S(0)[(i+w2-1)%w2+w2]);  // down-right

	 // nearest neighbor bonds, second column
	 H1 += inner(S(0)[i+w2], S(1)[i+w2]);         // horizontal
	 H1 += inner(S(0)[i+w2], S(1)[(i+1)%w2]);     // up-right
	 H1 += inner(S(0)[i+w2], S(1)[i]);            // down-right

	 // next-nearest neighbor bonds, first column
	 // vertical
	 H2 += inner(S(0)[i], S(0)[(i+1)%w2]);

	 // 60 degree bonds
	 H2 += inner(S(0)[i], S(1)[i+w2]);
	 H2 += inner(S(0)[i], S(1)[(i+w2-1)%w2+w2]);

	 // next-nearest neighbor bonds, second column
	 // vertical
	 H2 += inner(S(0)[i+w2], S(0)[(i+1)%w2+w2]);

	 // 60 degree bonds
	 H2 += inner(S(0)[i+w2], S(2)[(i+1)%w2]);
	 H2 += inner(S(0)[i+w2], S(2)[i]);

	 H_Jcell += inner(S(0)[i], S(0)[i+w2]);
	 H_Jcell += inner(S(0)[(i+1)%w2], S(0)[i+w2]);
      }


      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      Lattice["H_J1"]    = sum_unit(H1);
      Lattice["H_J2"]    = sum_unit(H2);
      Lattice["H_Jcell"] = sum_unit(H_Jcell);

      Lattice["Sa"] = sum_unit(S_A);
      Lattice["Sb"] = sum_unit(S_B);
      Lattice["Sc"] = sum_unit(S_C);
      Lattice["Stag_n60"] = sum_unit(S_stag_n60, w*2);

      Lattice.func("THM")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0)
         = "J1*H_J1 + J2*H_J2"; // an old lattice function, used in few projects in 2014-15.

      Lattice.func("THM2")(arg("J2") = 0.0)
              = "H_J1 + J2*H_J2";

      // Momentum operator in Y direction
      //Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);
      UnitCellMPO Ty = I(0);
      for (int i = 0; i < w-1; ++i)
      {
         Ty = Ty * Cell.swap_gate_no_sign((i/2)+(i%2)*w2, ((i+1)/2)+((i+1)%2)*w2);
      }
      Lattice["Ty"] = prod_unit_left_to_right(Ty.MPO(), w);

      // Reflection about X axis
      UnitCellMPO Ry = I(0);
      for (int i = 0; i < w2/2; ++i)
      {
	 Ry = Ry * Cell.swap_gate_no_sign(w2+i, w-i-1);
	 if (w2-i-1 > i+1)
	 {
	    Ry = Ry * Cell.swap_gate_no_sign(i+1, w2-i-1);
	 }
      }
      Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w);

      // for even size unit cell, add rotation by pi
      if (w2%2 == 0)
      {
	 UnitCellMPO TyPi = I(0);
	 for (int i = 0; i < w2/2; ++i)
	 {
	    TyPi = TyPi * Cell.swap_gate_no_sign(i, i+w2/2);
	    TyPi = TyPi * Cell.swap_gate_no_sign(i+w2, i+w2+w2/2);
	 }
	 Lattice["TyPi"] = prod_unit_left_to_right(TyPi.MPO(), w);
      }

      // Information about the lattice
      Lattice.set_description("SU(2) triangular lattice Heisenberg model, XC configuration");
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

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
