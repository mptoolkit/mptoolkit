// -*- C++ -*- $Id: bosehubbard-u1.cpp 1495 2015-05-20 16:11:06Z ianmcc $

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
      int MaxN = DefaultMaxN;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("NumBosons,N", prog_opt::value(&MaxN), 
	  FormatDefault("Maximum number of bosons per site", MaxN).c_str())
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
		   << "H_J    - nearest neighbor hopping\n"
		   << "H_U    - on-site Coulomb repulsion N*(N-1)/2\n"
	    //<< "\nOperator functions:\n"
	    //<< "H_flux{theta} - rung flux hopping cos(theta)*H_K + sin(theta)*H_Kc\n"
	    ;
         return 1;
      }

      LatticeSite Site = BosonU1(MaxN);
      UnitCell Cell(Site);
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2");

      InfiniteLattice Lattice(Cell);
      
      Lattice["H_J"] = sum_unit(BH(0)*B(1) + B(0)*BH(1));
      Lattice["H_U"] = sum_unit(0.5*N2(0));
      
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
