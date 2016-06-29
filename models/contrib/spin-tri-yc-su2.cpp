// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-tri-yc-su2.cpp
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

// Description: spin systems on triangular lattices with YC structure and efficient way of numbering; SU(2)-symmetric.
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
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>


namespace prog_opt = boost::program_options;


int IntPow(int x, int p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * IntPow(x, p-1);
}

// Functions that produces dimer terms for the Hamiltonian (including their h.c. (rotated) terms): 

UnitCellMPO DimerKinetic(half_int spin, int width, int cell1, int site1, int cell2, int site2, int cell3, int site3, int cell4, int site4) {

 LatticeSite Site = SpinSU2(spin);
 UnitCell Cell = repeat(Site, width);
 UnitCellOperator S(Cell, "S"), I(Cell, "I");
 UnitCellMPO H, Dots, DoubleDots, CrossDots;

 Dots = inner(S(cell1)[site1], S(cell4)[site4]) + inner(S(cell2)[site2], S(cell3)[site3]) - inner(S(cell1)[site1], S(cell2)[site2]) -   
        inner(S(cell3)[site3], S(cell4)[site4]) - inner(S(cell1)[site1], S(cell3)[site3]) - inner(S(cell2)[site2], S(cell4)[site4]);

 DoubleDots = inner(S(cell1)[site1], S(cell2)[site2])*inner(S(cell3)[site3], S(cell4)[site4]) +  
              inner(S(cell1)[site1], S(cell3)[site3])*inner(S(cell2)[site2], S(cell4)[site4]);

 CrossDots = inner(cross(S(cell1)[site1], S(cell2)[site2]), cross(S(cell3)[site3], S(cell4)[site4])) +
             inner(cross(S(cell1)[site1], S(cell3)[site3]), cross(S(cell2)[site2], S(cell4)[site4]));

 H = 0.125*I + 0.25*Dots + 0.5*DoubleDots - 0.125*CrossDots;

 return H;

}

UnitCellMPO DimerPotential(half_int spin, int width, int cell1, int site1, int cell2, int site2, int cell3, int site3, int cell4, int site4) {

 LatticeSite Site = SpinSU2(spin);
 UnitCell Cell = repeat(Site, width);
 UnitCellOperator S(Cell, "S"), I(Cell, "I");
 UnitCellMPO H, P12, P34, P13, P24;

 P12 = 0.25*I - inner(S(cell1)[site1], S(cell2)[site2]);
 P34 = 0.25*I - inner(S(cell3)[site3], S(cell4)[site4]);

 P13 = 0.25*I - inner(S(cell1)[site1], S(cell3)[site3]);
 P24 = 0.25*I - inner(S(cell2)[site2], S(cell4)[site4]);

 H = P12*P34 + P13*P24;

 return H;

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
	 std::cerr << "Constructs a triangular lattice in the YC configuration with wrapping vector (0,1) (in unit of\n"
		   << "+/-60 degree principal directions), while spins sit on edges. This creates an efficient way of\n" 
                   << "numbering for the 1D chain.\n\n"
                   << "Operators:\n"
		   << "H_J1     - nearest neighbor spin exchange\n"
		   << "H_J2     - next-nearest neighbor spin exchange\n"
                   << "H_chi    - explicit chiral term over NN triangular plaquettes\n" 
                   << "H_t      - kinetic term of quantum dimer model's Hamiltonian on the triangular lattice\n"
                   << "H_v      - potential term of quantum dimer model's Hamiltonian on the triangular lattice\n"
                   << "S        - total spin on a leg of the cylinder\n"
                   << "StagS    - staggered magnetization over a unit-cell\n"
                   << "Trans    - translation by one site (rotation by 2\u0071/w) in lattice short direction\n"
                   << "Ty       - momentum operator in lattice short direction\n"
                   << "Ref      - pure reflection in lattice short direction (may need applying T-operators to become" 
                   << "           general reflection)\n"
                   << "SwapWrap - changing the wraaping vector of lattice between 'old' and 'new' way of numbering"
		   //<< "*If the lattice could be potentially tripartite (width is a multiple of 3), then we\n"
		   //<< "define sublattice spin operators on a \"width*3\" unit cells as,\n"
		   << "S_A      - tripartite sublattice spin, including site S(0)[0]\n"
		   << "S_B      - tripartite sublattice spin, including site S(0)[1]\n"
		   << "S_C      - tripartite sublattice spin, including site S(0)[2]\n\n"
                   << "Functions:\n"
                   << "H2(J2 = NNN coupling strength, J_chi = chiral term coupling strength)\n\n"                         
	    ;
         return 1;
      }

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, w);
      UnitCellOperator S(Cell, "S"), StagS(Cell, "StagS");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      // Add some operators on the unit-cell
   
      for (int i = 0; i < w; ++i)
      {
	 S += S[i];                    // total spin on a leg of cylinder
         StagS += IntPow(-1,i) * S[i]; // only one of three possible staggered magnetization formation.            
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
       {
           //T *= 0.5*( 0.25*inner(S[i],S[i+1]) + 1 );
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
       }

      Ref = I(0);  // old way of representing an explicit R-operator.
      for (int i = 0; i < w/2; ++i)
       {
           //R *= 0.5*( 0.25*inner(S[i],S[w-i-1]) + 1 );
           Ref = Ref(0) * Cell.swap_gate_no_sign(i, w-i-1);
       }
          
      // if we could have tripartite symmetry, add operators for the sublattice magnetization
      UnitCellMPO S_A, S_B, S_C;

      //if (w%3 == 0)
      //{
	 for (int i = 0; i < w; i += 3)
	 {
	    S_A += S(0)[i]   + S(1)[(i+2)%w] + S(2)[(i+1)%w];
	    S_B += S(0)[i+1] + S(1)[(i+3)%w] + S(2)[(i+2)%w];
	    S_C += S(0)[i+2] + S(1)[(i+4)%w] + S(2)[(i+3)%w];
	 }
      //}


      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2, Ht, Hv, Hy, Hchi;
    
      for (int i = 0; i < w; ++i)
      {
	// nearest neighbor bonds
	 // vertical bonds
	 H1 += inner(S(0)[i], S(0)[(i+1)%w]);
	 // 60 degree bonds
	 H1 += inner(S(0)[i], S(1)[i]);
	 H1 += inner(S(0)[i], S(1)[(i+1)%w]);

	 // vertical bonds only
	 Hy += inner(S(0)[i], S(0)[(i+1)%w]);

	// next-nearest neighbor bonds
	 H2 += inner(S(0)[i], S(2)[(i+1)%w]);                  // horizontal
	 H2 += inner(S(0)[i], S(1)[(i+w-1)%w]);                // up-right
	 H2 += inner(S(0)[i], S(1)[(i+2)%w]);                  // down-right
         
        // kinetic terms of the dimer model
         Ht += DimerKinetic(Spin, w, 0, i, 1, i, 1, (i+1)%w, 2, (i+1)%w);               // horizontal rhombus terms
         Ht += DimerKinetic(Spin, w, 0, i, 0, (i+1)%w, 1, i, 1, (i+1)%w);               // upper-vertical rhombus terms
         Ht += DimerKinetic(Spin, w, 0, i, 1, (i+1)%w, 0, (i+1)%w, 1, (i+2)%w);         // lower-vertical rhombus terms

        // potential terms of the dimer model
         Hv += DimerPotential(Spin, w, 0, i, 1, i, 1, (i+1)%w, 2, (i+1)%w);             // horizontal rhombus terms
         Hv += DimerPotential(Spin, w, 0, i, 0, (i+1)%w, 1, i, 1, (i+1)%w);             // upper-vertical rhombus terms
         Hv += DimerPotential(Spin, w, 0, i, 1, (i+1)%w, 0, (i+1)%w, 1, (i+2)%w);       // lower-vertical rhombus terms

	 // right-facing triangle
	 Hchi += inner(S(0)[i], cross(S(0)[(i+1)%w], S(1)[(i+1)%w]));

	 // left-facing triangle
	 Hchi += inner(S(0)[i], cross(S(1)[(i+1)%w], S(1)[i]));
      }

      // Reflection.  This is in the 'wrong' 45 degree angle
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

      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      Lattice["H_J1"]  = sum_unit(H1);
      Lattice["H_y"]   = sum_unit(Hy);
      Lattice["H_J2"]  = sum_unit(H2);
      Lattice["H_t"]   = sum_unit(Ht);
      Lattice["H_v"]   = sum_unit(Hv);
      Lattice["H_chi"] = sum_unit(Hchi);

      Lattice.func("H")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0)
             = "J1*H_J1 + J2*H_J2"; // old lattice function 

      Lattice.func("H2")(arg("J2") = 0.0, arg("J_chi") = 0.0)
              = "H_J1 + J2*H_J2 + J_chi*H_chi";

      // Add the tripartite sublattice magnetization operators
      //if (w%3 == 0)
      //{
	 Lattice["Sa"] = sum_unit(S_A, w*3);
	 Lattice["Sb"] = sum_unit(S_B, w*3);
	 Lattice["Sc"] = sum_unit(S_C, w*3);
      //}

      // Momentum operators in Y-direction
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);
      Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w*w);

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

      // Reflection.  Fixed to reflect about a horizontal axis.  This is the reverse order of unit cells to the R-"45 degree" operator.
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

      // SwapWrap. Change between wrapping vectors of 'old' and 'new' way of numbering of the lattice.
      /* UnitCellMPO SwapWrap = I(0);
       for (int c = 0; c < w; ++c)
      {
         UnitCellMPO ThisWrap = I(0);

         if (c != 0)
         { 
           for (int i1 = 0; i1 < w-c; ++i1)
           {
              int i2 = (i1+w-c)%w;
              ThisWrap = ThisWrap * Cell.swap_gate_no_sign(i1,i2);
           }
         }

         ThisWrap.translate(c*w);
         SwapWrap = SwapWrap * ThisWrap;          
      } */         
      Lattice["SwapWrap"] = prod_unit_left_to_right(UnitCellMPO(Ref(0)).MPO(), w);  

      // 'identity' operator in the spin-1/2 auxiliary basis
      Lattice["I_2"] = prod_unit_left_to_right(UnitCellMPO(I(0)).MPO(), w)
	 * ProductMPO::make_identity(UnitCellMPO(I(0)).MPO().LocalBasis2List(), 
				     QuantumNumber(Cell.GetSymmetryList(), "0.5"));

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
