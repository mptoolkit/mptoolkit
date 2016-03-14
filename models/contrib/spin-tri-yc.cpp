// -*- C++ -*- $Id: spin-tri-yc.cpp 1490 2015-05-19 09:15:06Z seyed $
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
#include "models/spin.h"
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
                   << "Sx      - total spin on a leg of the cylinder\n"
                   << "Sy      - total spin on a leg of the cylinder\n"
                   << "Sz      - total spin on a leg of the cylinder\n"
                   << "Trans   - translation by one site (rotation by 2\u0071/w) in lattice short direction\n"
                   << "Ref     - reflection in lattice short direction (may need applying T-operators to become gneral reflection)\n"
                   << "Functions:\n"
                   << "not yet implemented ...\n\n"                          
	    ;
         return 1;
      }

      LatticeSite Site = SpinSite(Spin);
      UnitCell Cell = repeat(Site, w);
      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Trans(Cell, "Trans"), Ref(Cell, "Ref");
      UnitCellOperator RyUnit(Cell, "RyUnit");

      // Add some operators on the unit-cell
      
      for (int i = 0; i < w; ++i)
      {
         Sx += Sx[i];                                         // total spin (x-component) on a leg of cylinder
         Sy += Sy[i];                                         // total spin (y-component) on a leg of cylinder
         Sz += Sz[i];                                         // total spin (z-component) on a leg of cylinder

      }

      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
       {
           Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
       }

      Ref = I(0);
      for (int i = 0; i < w/2; ++i)
       {
           Ref = Ref(0) * Cell.swap_gate_no_sign(i, w-i-1);
       }
          
      // if we could have tripartite symmetry, add operators for the sublattice magnetization
      //UnitCellMPO S_A, S_B, S_C;


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

      // Construct the Hamiltonian for a single unit-cell,
      //UnitCellMPO H1, H2, Ht, Hv;

      //Lattice["H_J1"] = sum_unit(H1);

      //Lattice.func("H")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0) 
      // = "J1*H_J1 + J2*H_J2";


      // Momentum operator in Y direction
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
