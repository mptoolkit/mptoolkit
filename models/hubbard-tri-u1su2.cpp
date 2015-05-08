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
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
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
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Operators:\n"
		   << "H_t     - nearest-neighbor hopping between the apex sites of clusters\n"
		   << "H_t2    - next-nearest-neighbor hopping between the apex sites of clusters\n"
		   << "H_tc    - hopping inside the cluster\n"
		   << "H_tp    - nearest-neighbor hopping between non-apex sites of clusters\n"
		   << "H_U     - on-site Coulomb interaction n_up*n_down\n"
		   << "H_Us    - on-site Coulomb interaction (n_up-1/2)(n_down-1/2)\n"
	    ;
         return 1;
      }

      LatticeSite Site = FermionU1SU2();
      UnitCell Cell(repeat(Site, 3));
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), Pdouble(Cell, "Pdouble"),
	 Hu(Cell, "Hu"), N(Cell, "N"), R(Cell, "R"), S(Cell, "S");
      // parity operators
      // **NOTE** Currently the swap_gate doesn't do fermion signs, so we
      // need to add this by hand
      Cell["Parity"] = Cell.swap_gate_no_sign(0,2) * exp(math_const::pi*std::complex<double>(0,1)
							 *((N[0]+N[2])*N[1]+N[0]*N[2]));
      // Reflection operator.
      R = R[0]*R[1]*R[2];

      // some operators per unit cell
      N = N[0] + N[1] + N[2];
      S = S[0] + S[1] + S[2];

      InfiniteLattice Lattice(Cell);

      Lattice["H_t"]  = sum_unit(-(dot(CH(0)[1], C(1)[1]) + dot(C(0)[1], CH(1)[1])));
      Lattice["H_t2"] = sum_unit(-(dot(CH(0)[1], C(2)[1]) + dot(C(0)[1], CH(2)[1])));
      Lattice["H_tc"] = sum_unit(-(dot(CH(0)[0], C(0)[1]) + dot(C(0)[0], CH(0)[1])
				   + dot(CH(0)[1], C(0)[2]) + dot(C(0)[1], CH(0)[2])
				   + dot(CH(0)[0], C(0)[2]) + dot(C(0)[0], CH(0)[2])));
      Lattice["H_tp"] = sum_unit(-(dot(CH(0)[0], C(1)[0]) + dot(C(0)[0], CH(1)[0])
				   + dot(CH(0)[2], C(1)[2]) + dot(C(0)[2], CH(1)[2])));
      Lattice["H_U"]  = sum_unit(Pdouble(0)[0] + Pdouble(0)[1] + Pdouble(0)[2]);
      Lattice["H_Us"] = sum_unit(Hu(0)[0] + Hu(0)[1] + Hu(0)[2]);

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
