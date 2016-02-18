// -*- C++ -*- $Id: spin-tri-yc-u1.cpp 1490 2015-05-19 09:15:06Z ianmcc $
// Authors: Ian P. McCulloch and Seyed N. Saadatmand
// Contact: s.saadatmand@uq.edu.au
// <OBELIX> @ /data5/uqssaada/git/mptoolkit/models/contrib/spin-tri-yc-u1.cpp

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

      OperatorDescriptions OpDescriptions;
      OpDescriptions.add_operators()
         ("H_J1"                   , "nearest neighbor spin exchange")
         ("H_J2"                   , "next-nearest neighbor spin exchange")
         ("H_LongRangeIsing_inter" , "inter-cell interactions of long-range Ising model")
         ;

      OpDescriptions.add_functions()
         ("HS"                 , "Haldane-Shastry Hamiltonian with Sz*Sz interactions, parametized by 'lambda' (exponential decay as exp(-lambda*r))")
         ("LongRangeIsing_YC4" , "long-range Ising model on a 4-leg YC structure, parametized by 'lambda00', 'lambda01', and 'lambda02'")
         ;
      
      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Constructs a triangular lattice in the YC configuration with lattice vector (0,1).\n";
	 std::cerr << OpDescriptions << '\n';
	 std::cerr << "only if the lattice could be potentially tripartite (width is a multiple of 3):\n"
		   << "define sublattice spin operators on 'width*3' unit-cells.\n";
	    ;
         return 1;
      }

      const double PI = 4.0*std::atan(1.0);

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell = repeat(Site, w);

      std::cout << "Building all Hamiltonian operators:\n";

      // Add some operators on the unit cell
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I");

      for (int i = 0; i < w; ++i)
      {
	 Sp += Sp[i];     // total S+ on a leg of cylinder
	 Sm += Sm[i];     // total S- on a leg of cylinder
	 Sz += Sz[i];     // total Sz on a leg of cylinder
         std::cout << "... " << std::flush;
      }

      // Now we construct the InfiniteLattice,
      InfiniteLattice Lattice(Cell);

      // Construct the Hamiltonian for a single unit-cell,
      UnitCellMPO H1, H2, H_inter;
      for (int i = 0; i < w; ++i)
      {
	 // THM - nearest neighbor bonds
	 // --> vertical bonds
	 H1 += Sz(0)[i]*Sz(0)[(i+1)%w]
	    + 0.5 * (Sp(0)[i]*Sm(0)[(i+1)%w] + Sm(0)[i]*Sp(0)[(i+1)%w]);
         std::cout << ". " << std::flush;
	 // --> 60 degree bonds
	 H1 += Sz(0)[i]*Sz(1)[i] + 0.5 * (Sp(0)[i]*Sm(1)[i] + Sm(0)[i]*Sp(1)[i]);
	 H1 += Sz(0)[i]*Sz(1)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(1)[(i+1)%w] + Sm(0)[i]*Sp(1)[(i+1)%w]);
         std::cout << ".. " << std::flush; 

	 // THM - next-nearest neighbor bonds
	 H2 += Sz(0)[i]*Sz(2)[(i+1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
	 H2 += Sz(0)[i]*Sz(1)[(i+w-1)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
	 H2 += Sz(0)[i]*Sz(1)[(i+2)%w] + 0.5 * (Sp(0)[i]*Sm(2)[(i+1)%w] + Sm(0)[i]*Sp(2)[(i+1)%w]);
         std::cout << "... " << std::flush;     

        for (int j = i+1; j < w; ++j)
        {
          // Long-range Ising - inter-cell interations
          H_inter += ( std::sin( (j-i)*PI/w ) / std::sin( PI/w ) ) * Sz(0)[i]*Sz(0)[j];
          std::cout << ". " << std::flush;
        }
      }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);
      Lattice["H_LongRangeIsing_inter"] = sum_unit(H_inter);
      std::cout << "... " << std::flush;

      Lattice.func("THM")(arg("J1") = "cos(theta)", arg("J2") = "sin(theta)", arg("theta") = "atan(alpha)", arg("alpha") = 0.0)
                  = "J1*H_J1 + J2*H_J2";
      std::cout << ". " << std::flush;

      // a basic function for Haldane-Shastry model with Sz*Sz interations
      Lattice.func("HS")(arg("lambda") = 0.5, arg("i") = "0", arg("j") = "0")
                  = "exp(-lambda)*sum_string_inner( S(0)[i], exp(-lambda)*I(0), S(0)[j] )";
      std::cout << ". " << std::flush;

      Lattice.func("LongRangeIsing_YC4")(arg("lambda00") = 0.5, arg("lambda01") = 0.5, arg("lambda02") = 0.5)
                  = "HS{lambda00,0,0} + HS{lambda00,1,1} + HS{lambda00,2,2} + HS{lambda00,3,3} + HS{lambda01,0,1} + HS{lambda01,1,2} + HS{lambda01,2,3} + HS{lambda01,3,1} + HS{lambda02,0,2} + HS{lambda02,1,3} + H_LongRangeIsing_inter";
      std::cout << ". " << std::flush;

      std::cout << ">>>finished.\n";

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
