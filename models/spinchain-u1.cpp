// -*- C++ -*- $Id$

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
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Operators:\n"
		   << "H_J1z   - nearest neighbor spin coupling Sz Sz\n"
		   << "H_J1t   - nearest neighbor spin exchange Sp Sm + Sm Sp\n"
		   << "H_J1    - nearest neighbor spin exchange = H_J1z + H_J1t\n"
		   << "H_J2z   - next-nearest neighbor spin coupling Sz Sz\n"
		   << "H_J2t   - next-nearest neighbor spin exchange Sp Sm + Sm Sp\n"
		   << "H_J2    - next-nearest neighbor spin exchange = H_J1z + H_J1t\n"
		   << "H_B1    - nearest neighbor biquadratic spin exchange (S.S)^2\n"
		   << "H_B2    - next-nearest neighbor biquadratic spin exchange (S.S)^2\n"
		   << "H_mu    - single-ion anistotropy, H_mu = sum_i Sz(i)^2\n"
	    ;
         return 1;
      }

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(Site);
      InfiniteLattice Lattice(Cell);
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");

      Lattice["H_J1z"] = sum_unit(Sz(0)*Sz(1));
      Lattice["H_J1t"] = sum_unit(0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)));
      Lattice["H_J1"] = sum_unit(Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)));

      Lattice["H_J2z"] = sum_unit(Sz(0)*Sz(2));
      Lattice["H_J2t"] = sum_unit(0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));
      Lattice["H_J2"] = sum_unit(Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));

      Lattice["H_B1"] = sum_unit(pow(Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)), 2));
      Lattice["H_B2"] = sum_unit(pow(Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)), 2));

      Lattice["H_mu"] = sum_unit(Sz(0)*Sz(0));

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
