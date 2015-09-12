// -*- C++ -*- $Id: spinchain-u1.cpp 1595 2015-09-04 13:15:02Z ianmcc $

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
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
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
	 ("H_J1z", "nearest neighbor spin coupling Sz Sz")
	 ("H_J1t", "nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
	 ("H_J1" , "nearest neighbor spin exchange = H_J1z + H_J1t")
	 ("H_J2z", "next-nearest neighbor spin coupling Sz Sz")
	 ("H_J2t", "next-nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
	 ("H_J2" , "next-nearest neighbor spin exchange = H_J1z + H_J1t")
	 ("H_B1" , "nearest neighbor biquadratic spin exchange (S.S)^2")
	 ("H_B2" , "next-nearest neighbor biquadratic spin exchange (S.S)^2")
	 ("H_mu" , "single-ion anistotropy, H_mu = sum_i Sz(i)^2")
	 ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Operators:\n" << OpDescriptions;
	 std::cerr << "only for spin-1: H_AKLT  - AKLT Hamiltonian H+J1 + (1/3)*H_B1\n";
         return 1;
      }

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(Site);
      InfiniteLattice Lattice(Cell);
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");

      Lattice["H_J1z"] = sum_unit(Sz(0)*Sz(1));
      Lattice["H_J1t"] = sum_unit((sqrt(0.5)*Sp(0))*(sqrt(0.5)*Sm(1)) + (sqrt(0.5)*Sm(0))*(sqrt(0.5)*Sp(1)));
      Lattice["H_J1"]  = Lattice["H_J1z"] + Lattice["H_J1t"];
      //sum_unit(Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)));

      Lattice["H_J2z"] = sum_unit(Sz(0)*Sz(2));
      Lattice["H_J2t"] = sum_unit(0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));
      Lattice["H_J2"] = sum_unit(Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));

      Lattice["H_B1"] = sum_unit(pow(Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)), 2));
      Lattice["H_B2"] = sum_unit(pow(Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)), 2));

      Lattice["H_mu"] = sum_unit(Sz(0)*Sz(0));

      if (Spin == 1)
      {
	 Lattice["H_AKLT"] = Lattice["H_J1"] + (1.0/3.0)*Lattice["H_B1"];
	 Lattice["H_AKLT"].set_description("AKLT Hamiltonian H_J1 + (1/3)*H_B1");
      }

      // Information about the lattice
      Lattice.set_description("U(1) Spin chain");
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice to disk
      pheap::ExportObject(FileName, Lattice);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
