// -*- C++ -*- $Id$

//
// YC configuration of a triangular lattice
// The unit cell is size Width*2.
//
// Example for a width-5 lattice (site numbers in brackets are periodic repeats in the vertical
// direction (ie top-left (0) is the same site as the bottom-left 0).
// Sites 5,6,7,8,9 are the second unit cell, eg 5 is (1)[0]
//
//    (0)
//  (5)|
//(0)| 4
// | 9<|
// 4<| 3
// | 8<|
// 3<| 2
// | 7<|
// 2<| 1
// | 6<|
// 1<| 0
// | 5<|
// 0<|(4)
// |(9)
//(4)
//

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
      UnitCellOperator S(Cell, "S");

      // Add some operators on the unit cell
      for (int i = 0; i < w; ++i)
      {
	 S += S[i];     // total spin
      }

      // if we have tripartite symmetry, add operators for the sublattice magnetization
      UnitCellMPO S_A, S_B, S_C;
      if (w%3 == 0)
      {
	 for (int i = 0; i < w; i += 3)
	 {
	    S_A += S(0)[i]   + S(1)[(i+1)%w] + S(2)[(i+2)%w];
	    S_B += S(0)[i+1] + S(1)[(i+2)%w] + S(2)[i%w];
	    S_C += S(0)[i+2] + S(1)[i%w]     + S(2)[(i+1)%w];
	 }
      }

      // Now we construct the InfiniteLattice
      InfiniteLattice Lattice(Cell);

      // Construct the Hamiltonian for a single unit cell
      UnitCellMPO H1, H2;
      for (int i = 0; i < w; ++i)
      {
	 // Nearest neighbor bonds
	 // vertical bonds
	 H1 += inner(S(0)[i], S(0)[(i+1)%w]);
	 // 60 degree bonds
	 H1 += inner(S(0)[i], S(1)[i]);
	 H1 += inner(S(0)[i], S(1)[(i+w-1)%w]);

	 // next-nearest neighbor bonds
	 H2 += inner(S(0)[i], S(2)[i]);             // horizontal
	 H2 += inner(S(0)[i], S(1)[(i+1)%w]);       // up-right
	 H2 += inner(S(0)[i], S(1)[(i+w-2)%w]);     // down-right   
      }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);

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
