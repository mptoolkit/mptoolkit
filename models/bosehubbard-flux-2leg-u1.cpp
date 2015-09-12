// -*- C++ -*- $Id: bosehubbard-flux-u1.cpp 1477 2015-05-12 12:56:57Z ianmcc $

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/boson-u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include <boost/program_options.hpp>

//
// b^\dagger_{0,r} = exp[-ir\phi] a^\dagger_{0,r}
// b^\dagger_{1,r} = a^\dagger_{1,r}
//
// Hamiltonian:
// -J \sum_{r} ( exp[i\phi] b^\dagger_{0,r+1} b_{0,r}
//                  + b^\dagger_{1,r+1} b_{1,r} + H.c )
// -J_\perp \sum{r} b^\dagger_{0,r} b_{1,r} + H.c.
// + U/2 \sum_{l,r} n_{l,r} (n_{l,r}-1)
//
// For parameters J, K, U, alpha, the Hamiltonian is:
// J*(H_J1 + cos(pi*alpha)*H_J0 + sin(pi*alpha)*H_Jc0) + K*H_K + U*H_U
//
// Current operator is derivative with respect to flux, so
//
// j_c = -sin(pi*alpha)*H_J0 + cos(pi*alpha)*H_Jc0
//
// j^\parallel_0,r = i * (e^{i\phi} a^\dagger_{0,r+1} a_{0,r} - e^{-i\phi} a^\dagger_{0,r} a_{0,r+1})
// j^\parallel_1,r = i * (a^\dagger_{0,r+1} a_{0,r} - a^\dagger_{0,r} a_{0,r+1})
//
// 
// \sum_r j^\parallel_0,r = i \cos \phi [ a^\dagger_{0,r+1} a_{0,r} - a^\dagger_{0,r} a_{0,r+1} ]
//                            - \sin \phi [  a^\dagger_{0,r+1} a_{0,r} + a^\dagger_{0,r} a_{0,r+1} ]
//                        = \cos \phi H_Jc0 - \sin \phi H_J0
// \sum_r j^\parallel_1,r = H_Jc1


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
      
      OperatorDescriptions OpDescriptions;
      OpDescriptions.add_operators()
	 ("H_J0",  "nearest neighbor leg 0 hopping")
	 ("H_J1",  "nearest neighbor leg 1 hopping")
	 ("H_J",   "nearest neighbor leg hopping = H_J0 + H_J1")
	 ("H_Jc0", "nearest neighbor leg 0 current")
	 ("H_Jc1", "nearest neighbor leg 1 current")
	 ("H_Jc",  "edge current H_Jc0 - H_Jc1")
	 ("H_K",   "nearest neighbor rung hopping")
	 ("H_Kc",  "nearest neighbor rung current")
	 ("H_U",   "on-site Coulomb repulsion N*(N-1)/2")
	 ;

      OpDescriptions.add_functions()
	 ("H", "Hamiltonian, parametized by K, alpha (flux), U, J")
	 ("jTotal", "Total current, function of alpha")
	 ("jTotalU0", "Current at U=0 as a function of K, alpha (number, not an operator)")
	 ;


      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << OpDescriptions;
         return 1;
      }

      LatticeSite Site = BosonU1(MaxN);
      UnitCell Cell(repeat(Site, 2));
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2"),
	 j0(Cell, "j0"), j1(Cell, "j1");

      j0 = std::complex<double>(0,1)*(BH(0)[0]*B(1)[0] - B(0)[0]*BH(1)[0]);
      j1 = std::complex<double>(0,1)*(BH(0)[1]*B(1)[1] - B(0)[1]*BH(1)[1]);

      Cell.func("jk")("K", "alpha") = "K*(BH(0)[0]*B(0)[1] - B(0)[0]*BH(0)[1]";
      Cell.func("Jp0")("K", "alpha") = "-cos(pi*alpha)*(BH(0)[0]*B(1)[0] + B(0)[0]*BH(1)[0])"
	 " + i*sin(pi*alpha)*(BH(0)[0]*B(1)[0] - B(0)[0]*BH(1)[0])";
      Cell.func("Jp1")("K", "alpha") = "-cos(pi*alpha)*(BH(0)[1]*B(1)[1] + B(0)[1]*BH(1)[1])"
	 " + i*sin(pi*alpha)*(BH(0)[1]*B(1)[1] - B(0)[1]*BH(1)[1])";

      Cell.func("jTotalU0")("K", "alpha", arg("J")=1.0) = "K^2 * sin(pi*alpha) / (8*J * "
      	 "sin(pi*alpha/2)^4 * (1+(K/(2*J*sin(pi*alpha/2)))^2)^0.5)";

      InfiniteLattice Lattice(Cell);
      
      UnitCellMPO HJ0, HJ1, HJc0, HJc1, HK, HKc, HU;
      HJ0 = BH(0)[0]*B(1)[0] + B(0)[0]*BH(1)[0];
      HJ1 = BH(0)[1]*B(1)[1] + B(0)[1]*BH(1)[1];
      HJc0 = std::complex<double>(0,1)*(BH(0)[0]*B(1)[0] - B(0)[0]*BH(1)[0]);
      HJc1 = std::complex<double>(0,1)*(BH(0)[1]*B(1)[1] - B(0)[1]*BH(1)[1]);
      HU = 0.5*(N2(0)[0] + N2(0)[1]);

      HK = BH(0)[0]*B(0)[1] + B(0)[0]*BH(0)[1];
      HKc = std::complex<double>(0,1)*(BH(0)[0]*B(0)[1] - B(0)[0]*BH(0)[1]);

      Lattice["H_J0"] = sum_unit(HJ0);
      Lattice["H_J1"] = sum_unit(HJ1);
      Lattice["H_J"] = sum_unit(HJ0 + HJ1);
      Lattice["H_Jc0"] = sum_unit(j0);
      Lattice["H_Jc1"] = sum_unit(j1);
      Lattice["H_Jc"] = sum_unit(j0-j1);
      Lattice["H_K"] = sum_unit(HK);
      Lattice["H_Kc"] = sum_unit(HKc);
      Lattice["H_U"] = sum_unit(HU);
      
      // For parameters J, K = J_\perp, U, alpha, the Hamiltonian is:
      // J*(H_J1 + cos(pi*alpha)*(H_J2 + H_J0) + sin(pi*alpha)*(H_Jc2 - H_Jc0)) + K*H_K + U*H_U

      Lattice.func("H")("K", arg("alpha")=0.0, arg("U")=0.0, arg("J")=1.0)
	 = "J*(H_J1 + cos(pi*alpha)*H_J0 + sin(pi*alpha)*H_Jc0) + K*H_K + U*H_U";

      Lattice.func("jTotal")("alpha") = "-sin(pi*alpha)*H_J0 + cos(pi*alpha)*H_Jc0";

      Lattice.set_description("Bosonic 2-leg ladder with flux");
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
