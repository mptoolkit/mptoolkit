// -*- C++ -*- $Id: spin-tri-yc-su2.EfficientNumbering.cpp 1490 2015-05-19 09:15:06Z ianmcc $
// Authors: Ian P. McCulloch and Seyed N. Saadatmand
// Contact: s.saadatmand@uq.edu.au

//
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
// |  8< |
// 2< | 15
// |  9< |
// 3< | 16
// | 10< |
// 4< | 17
// | 11< |
// 5< | (12)
// | (6)
//(0)
//

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-u1.h"
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
	 std::cerr << "Constructs a triangular lattice in the YC configuration with lattice vector (0,1) and\n"
		   << "efficient way of numbering 1D chain. Operators are:\n"
		   << "H_J1    - nearest neighbor spin exchange\n"
		   << "H_J2    - next-nearest neighbor spin exchange\n"
		   << "\nIf the lattice could potentially be tripartite (width is a multiple of 3), then we\n"
		   << "define sublattice spin operators on width*3 unit cells,\n"

	    ;
         return 1;
      }

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell = repeat(Site, w);

      // Add some operators on the unit cell
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      for (int i = 0; i < w; ++i)
      {
	 Sp += Sp[i];     // total spin on a leg of cylinder
	 Sm += Sm[i];     // total spin on a leg of cylinder
	 Sz += Sz[i];     // total spin on a leg of cylinder
      }

      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2;
      for (int i = 0; i < w; ++i)
      {
	 // Nearest neighbor bonds
	 // vertical bonds
	 H1 += Sz(0)[i]*Sz(0)[(i+1)%w]
	    + 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%w] + Sm(0)[i]*Sp(0)[(i+1)%w]);
	 // 60 degree bonds
	 H1 += Sz(0)[i]*Sz(1)[i] + 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
	 H1 += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+1)%w] + Sm(0)[i]*Sp(1)[(i+1)%w]);

	 // next-nearest neighbor bonds
	 H2 += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
	 H2 += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
	 H2 += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
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
