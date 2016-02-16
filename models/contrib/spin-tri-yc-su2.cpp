// -*- C++ -*-
// Description: spin systems on triangular lattices with YC structure and efficient way of numbering; SU(2)-symmetric. spin-tri-yc-su2.cpp   
// Authors: Ian P. McCulloch and Seyed N. Saadatmand
// Contact: s.saadatmand@uq.edu.au

// <obelix> @ /data5/uqssaada/git/mptoolkit/models/contrib/

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

// Functions that produces dimer terms for the Hamiltonian (including their h.c. or rotated term): 

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
	 std::cerr << "Constructs a triangular lattice in the YC configuration with wrapping vector (0,1),\n"
		   << "while spins sit on edges. This employs an efficient way of numbering in the 1D chain.\n\n"
                   << "Operators:\n"
		   << "H_J1    - nearest neighbor spin exchange\n"
		   << "H_J2    - next-nearest neighbor spin exchange\n"
                   << "H_t     - kinetic term of quantum dimer model's Hamiltonian on the triangular lattice\n"
                   << "H_v     - potential term of quantum dimer model's Hamiltonian on the triangular lattice\n"
                   << "S       - total spin on a leg of the cylinder\n"
                   << "StagS   - staggered magnetization over a unit-cell\n\n"
                   << "Trans   - translation by one site (rotation by 2\u0071/w) in lattice short direction\n"
                   << "Ref     - reflection in lattice short direction (may need applying T-operators to become gneral reflection)\n"
		   << "If the lattice could be potentially tripartite (width is a multiple of 3), then we\n"
		   << "define sublattice spin operators on a \"width*3\" unit cells as,\n"
		   << "S_A     - tripartite sublattice spin, including site S(0)[0]\n"
		   << "S_B     - tripartite sublattice spin, including site S(0)[1]\n"
		   << "S_C     - tripartite sublattice spin, including site S(0)[2]\n\n"
                   << "Functions:\n"
                   << "H( J1 = NN coupling strength, J2 = NNN coupling strength, theta = atan(J2/J1)\n" 
                   << "  \"radians\", alpha = J2/J1 )\n\n"                         
	    ;
         return 1;
      }

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell = repeat(Site, w);
      UnitCellOperator S(Cell, "S"), StagS(Cell, "StagS");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");

      // Add some operators on the unit-cell
   
      for (int i = 0; i < w; ++i)
      {
	 S += S[i];                                         // total spin on a leg of cylinder
         StagS += IntPow(-1,i) * S[i];            
      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
       {
           //T *= 0.5*( 0.25*inner(S[i],S[i+1]) + 1 );
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
       }

      Ref = I(0);
      for (int i = 0; i < w/2; ++i)
       {
           //R *= 0.5*( 0.25*inner(S[i],S[w-i-1]) + 1 );    // old way of representing R-operator.
           Ref = Ref(0) * Cell.swap_gate_no_sign(i, w-i-1);
       }
          
      // if we could have tripartite symmetry, add operators for the sublattice magnetization
      UnitCellMPO S_A, S_B, S_C;

      if (w%3 == 0)
      {
	 for (int i = 0; i < w; i += 3)
	 {
	    S_A += S(0)[i]   + S(1)[(i+2)%w] + S(2)[(i+1)%w];
	    S_B += S(0)[i+1] + S(1)[(i+3)%w] + S(2)[(i+2)%w];
	    S_C += S(0)[i+2] + S(1)[(i+4)%w] + S(2)[(i+3)%w];
	 }
      }


      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2, Ht, Hv, Hy;
    
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
      }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_y"] = sum_unit(Hy);
      Lattice["H_J2"] = sum_unit(H2);
      Lattice["H_t"] = sum_unit(Ht);
      Lattice["H_v"] = sum_unit(Hv);

      Lattice.func("H")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0)
	 = "J1*H_J1 + J2*H_J2";

      // Add the tripartite sublattice magnetization operators
      if (w%3 == 0)
      {
	 Lattice["S_A"] = sum_unit(S_A, w*3);
	 Lattice["S_B"] = sum_unit(S_B, w*3);
	 Lattice["S_C"] = sum_unit(S_C, w*3);
      }

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
