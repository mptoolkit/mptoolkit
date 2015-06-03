// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/boson-u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      int MaxN = 5;
      int Width = 2;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("NumBosons,N", prog_opt::value(&MaxN), 
	  FormatDefault("Maximum number of bosons per site", MaxN).c_str())
	 ("width,w", prog_opt::value(&Width),
	  FormatDefault("Width of the ladder", Width).c_str())
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
	 std::cerr << "Operators:\n"
		   << "H_J   - nearest neighbor hopping on the legs\n"
		   << "H_K   - nearest neighbor hopping on the rungs\n"
		   << "H_U   - on-site Coulomb repulsion N*(N-1)/2\n"
		   << "H_U12 - nearest-neighbor Coulomb repulsion on the rungs\n"
	    ;
         return 1;
      }

      LatticeSite Site = BosonU1(MaxN);
      UnitCell Cell(repeat(Site, Width));
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2"),
	 Delta(Cell, "D");

      // difference of boson occupation numbers between edges of the ladder
      Delta = N(0)[Width-1] - N(0)[0];

      InfiniteLattice Lattice(Cell);
      
      UnitCellMPO HJ, HK, HU, HU12;
      for (int i = 0; i < Width; ++i)
      {
	 HJ -= BH(0)[i]*B(1)[i] + B(0)[i]*BH(1)[i];
      }

      for (int i = 0; i < Width-1; ++i)
      {
	 HK -= BH(0)[i]*B(0)[i+1] + B(0)[i]*BH(0)[i+1];
	 HU += 0.5*N2(0)[i];
      }
      HU12 = N(0)[0] * N(0)[Width-1];

      Lattice["H_J"]   = sum_unit(HJ);
      Lattice["H_K"]   = sum_unit(HK);
      Lattice["H_U"]   = sum_unit(HU);
      Lattice["H_U12"] = sum_unit(HU12);
      Lattice["D"]     = sum_unit(Delta);
      Lattice["D2"]    = sum_unit(Delta*Delta);
      
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
