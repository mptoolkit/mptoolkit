// -*- C++ -*- $Id$

//
// YC configuration of a triangular lattice
// The unit cell is size Width*2.
//
// Example for a width-5 lattice (site numbers in brackets are periodic repeats)
//
//(0) (0)
// |(5)|
// 4<| 4
// | 9<|
// 3<| 3
// | 8<|
// 2<| 2
// | 7<|
// 1<| 1
// | 6<|
// 0<| 0
// | 5<|
//(4)|(4)
//  (9)
//
// A B A
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
      int Width = 4;
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
	 ("width,w", prog_opt::value(&Width), "width of the cylinder [default 4]")
         ("out,o", prog_opt::value(&LatticeName), "output filename [required]")
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
         std::cerr << "usage: spin-yc-su2 [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Constructs a triangular lattice in the YC configuration with lattice vector (0,1)\n"
		   << "Operators:\n"
		   << "H_J1    - nearest neighbor spin exchange\n"
		   << "H_J2    - next-nearest neighbor spin exchange\n"
	    ;
         return 1;
      }

      LatticeSite Site = CreateSU2SpinSite(Spin);
      UnitCell Cell = repeat(Site, Width*2);
      UnitCellOperator S(Cell, "S");

      // Add some operators on the unit cell
      for (int i = 0; i < Width*2; ++i)
      {
	 S += S[i];     // total spin
      }

      InfiniteLattice Lattice(Cell);

      // Construct the Hamiltonian for a single unit cell
      UnitCellMPO H1, H2;
      for (int i = 0; i < Width; ++i)
      {
	 // Nearest neighbor bonds
	 // Strip A
	 // vertical bonds
	 H1 += inner(S(0)[i], S(0)[(i+1)%Width]);
	 // 60 degree bonds
	 H1 += inner(S(0)[i], S(0)[i+Width]);
	 H1 += inner(S(0)[i], S(0)[(i+1)%Width+Width]);

	 // strip B
	 // vertical bonds
	 H1 += inner(S(0)[i+Width], S(0)[(i+1)%Width+Width]);
	 // 60 degree bonds
	 H1 += inner(S(0)[i+Width], S(1)[i]);
	 H1 += inner(S(0)[i+Width], S(1)[(i+Width-1)%Width]);

	 // next-nearest neighbor bonds
	 H2 += inner(S(0)[i], S(1)[i]);
	 H2 += inner(S(0)[i], S(0)[(i+2)%Width+Width]);
	 H2 += inner(S(0)[i], S(0)[(i+Width-1)%Width+Width]);

	 H2 += inner(S(0)[i+Width], S(1)[i+Width]);
	 H2 += inner(S(0)[i+Width], S(1)[(i+1)%Width]);
	 H2 += inner(S(0)[i+Width], S(1)[(i+Width-2)%Width]);
      }

      Lattice["H_J1"] = sum_unit(H1);
      Lattice["H_J2"] = sum_unit(H2);

      TRACE(Lattice.GetUnitCell().size());

      pheap::ExportObject(LatticeName, Lattice);
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
