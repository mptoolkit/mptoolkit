// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int x = 0;
      int y = 4;
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the spin", Spin).c_str())
	 (",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str())
	 (",y", prog_opt::value(&y), FormatDefault("y wrapping vector", y).c_str())
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
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Spin cylinder.Wrapping vector is (x,y).  y is the width, x is the offset.  (0,y) is the\n"
		   << "YC configuration with width and unit cell size y.  (x,y) is with y cylinder\n"
		   << "with offset x.  The unit cell size is x.\n";
	 std::cerr << "Operators:\n"
		   << "H_J1x   - nearest neighbor spin exchange in the x direction\n"
		   << "H_J1y   - nearest neighbor spin exchange in the y direction\n"
		   << "H_J1    - nearest neighbor spin exchange H_J1x + H_J1y\n"
	    ;
         return 1;
      }

      int CellSize = x == 0 ? y : x;

      LatticeSite Site = SpinSU2(Spin);
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(Cell);
      UnitCellOperator S(Cell, "S"), Q(Cell, "Q");

      UnitCellMPO J1x, J1y;
      // the XY configuration is special
      if (x == 0)
      {
	 for (int i = 0; i < y; ++i)
	 {
	    J1x += inner(S(0)[i], S(1)[i]);
	    J1y += inner(S(0)[i], S(0)[(i+1)%y]);
	 }
      }
      else
      {
	 for (int i = 0; i < x-1; ++i)
	 {
	    J1x += inner(S(0)[i], S(0)[i+1]);
	 }
	 J1x += inner(S(0)[x-1], S(y+1)[0]);
	 for (int i = 0; i < x; ++i)
	 {
	    J1y += inner(S(0)[i], S(1)[i]);
	 }
      }

      Lattice["H_J1x"] = sum_unit(J1x);
      Lattice["H_J1y"] = sum_unit(J1y);
      Lattice["H_J1"] = sum_unit(J1x+J1y);

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
