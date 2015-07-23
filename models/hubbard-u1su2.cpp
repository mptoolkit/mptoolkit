// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/fermion-u1su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
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
	 ("H_t"  , "nearest neighbor hopping")
	 ("H_t2" , "next-nearest neighbor hopping")
	 ("H_U"  , "on-site Coulomb interaction n_up*n_down")
	 ("H_Us" , "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
	 ("H_V"  , "nearest-neighbor Coulomb interaction")
	 ("H_J"  , "nearest-neighbor complex hopping i*(C^\\dagger_i C_{i+1} - H.c.)")
	 ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      LatticeSite Site = FermionU1SU2();
      UnitCell Cell(Site);
      InfiniteLattice Lattice(Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), Pdouble(Cell, "Pdouble"),
	 Hu(Cell, "Hu"), N(Cell, "N");

      Lattice["H_t"]  = sum_unit(dot(CH(0), C(1)) + dot(C(0), CH(1)));
      Lattice["H_t2"] = sum_unit(dot(CH(0), C(2)) + dot(C(0), CH(2)));
      Lattice["H_U"]  = sum_unit(Pdouble(0));
      Lattice["H_Us"] = sum_unit(Hu(0));
      Lattice["H_V"]  = sum_unit(dot(N(0), N(1)));
      Lattice["H_J"]  = sum_unit(std::complex<double>(0,1)
				 *(dot(CH(0), C(1)) - dot(C(0), CH(1))));

      // Information about the lattice
      Lattice.set_description("U(1)xSU(2) Fermi Hubbard model");
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

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
