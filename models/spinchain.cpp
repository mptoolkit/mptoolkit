// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int S = 0.5;
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&S), "magnitude of the spin [default 0.5]")
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
		   << "H_xx    - nearest neighbor spin coupling Sx Sx\n"
		   << "H_yy    - nearest neighbor spin exchange Sy Sy\n"
		   << "H_zz    - nearest neighbor spin exchange Sz Sz\n"
		   << "H_x     - magnetic field in the x direction\n"
		   << "H_y     - magnetic field in the y direction\n"
		   << "H_z     - magnetic field in the z direction\n"

	    ;
         return 1;
      }

      LatticeSite Site = SpinSite(S);
      UnitCell Cell(Site);
      InfiniteLattice Lattice("Spin chain", Cell);
      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz");

      Lattice["H_xx"] = sum_unit(Sx(0)*Sx(1));
      Lattice["H_yy"] = sum_unit(Sy(0)*Sy(1));
      Lattice["H_zz"] = sum_unit(Sz(0)*Sz(1));
      Lattice["H_x"] = sum_unit(Sx(0));
      Lattice["H_y"] = sum_unit(Sy(0));
      Lattice["H_z"] = sum_unit(Sz(0));

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
