// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/boson-2component-u1z2.h"
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
		   << "H_J   - nearest neighbor hopping on the legs\n"
		   << "H_K   - tunneling between components\n"
		   << "H_U   - inter-species Coulomb repulsion\n"
		   << "H_U12 - intra-species Coulomb repulsion\n"
	    ;
         return 1;
      }

      LatticeSite Site = Boson2ComponentU1Z2(MaxN);
      UnitCell Cell = Site;
      UnitCellOperator BH_A(Cell, "BH_A"), B_A(Cell, "B_A"), N_A(Cell, "N_A"), N2_A(Cell, "N2_A"),
	 BH_S(Cell, "BH_S"), B_S(Cell, "B_S"), N_S(Cell, "N_S"), N2_S(Cell, "N2_S");

      InfiniteLattice Lattice(Cell);
      
      UnitCellMPO HJ = -(BH_A(0)*B_A(1) + B_A(0)*BH_A(1) + BH_S(0)*B_S(1) + B_S(0)*BH_S(1));
      UnitCellMPO HK = -(N_S(0) - N_A(0));
      
      UnitCellMPO PairHopping = pow(BH_S(0)*B_A(0),2) + pow(BH_A(0)*B_S(0),2);
      
      UnitCellMPO HU = N_S(0)*N_A(0) + 0.25 * (N2_S(0) + N2_A(0) + PairHopping);
      UnitCellMPO HU12 = 0.25 * (N2_S(0) + N2_A(0) - PairHopping);

      Lattice["H_J"]   = sum_unit(HJ);
      Lattice["H_K"]   = sum_unit(HK);
      Lattice["H_U"]   = sum_unit(HU);
      Lattice["H_U12"] = sum_unit(HU12);
      
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
