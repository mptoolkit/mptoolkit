// -*- C++ -*- $Id$

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
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
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
         std::cerr << "usage: spinchain-su2 [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Operators:\n"
		   << "H_J1    - nearest neighbor spin exchange\n"
		   << "H_J2    - next-nearest neighbor spin exchange\n"
		   << "H_J3    - next-next-nearest neighbor spin exchange\n"
		   << "H_B1    - nearest neighbor biquadratic spin exchange (S.S)^2\n"
		   << "H_B2    - next-nearest neighbor biquadratic spin exchange (S.S)^2\n"
		   << "H_B3    - next-next-nearest neighbor biquadratic spin exchange (S.S)^2\n"
		   << "H_Q1    - nearest neighbor quadrupole exchange (Q.Q)\n"
		   << "H_Q2    - next-nearest neighbor quadrupole exchange (Q.Q)\n"
		   << "H_Q3    - next-next-nearest neighbor quadrupole exchange (Q.Q)\n"
	    ;
         return 1;
      }

      LatticeSite Site = CreateSU2SpinSite(Spin);
      UnitCell Cell(Site);
      InfiniteLattice Lattice(Cell);
      UnitCellOperator S(Cell, "S"), Q(Cell, "Q");

      Lattice["H_J1"] = sum_unit(inner(S(0), S(1)));
      Lattice["H_J2"] = sum_unit(inner(S(0), S(2)));
      Lattice["H_J3"] = sum_unit(inner(S(0), S(3)));

      Lattice["H_B1"] = sum_unit(pow(inner(S(0), S(1)), 2));
      Lattice["H_J2"] = sum_unit(pow(inner(S(0), S(2)), 2));
      Lattice["H_J3"] = sum_unit(pow(inner(S(0), S(3)), 2));

      Lattice["H_Q1"] = sum_unit(inner(Q(0), Q(1)));
      Lattice["H_Q2"] = sum_unit(inner(Q(0), Q(2)));
      Lattice["H_Q3"] = sum_unit(inner(Q(0), Q(3)));

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
