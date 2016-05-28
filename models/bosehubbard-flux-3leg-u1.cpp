// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/bosehubbard-flux-3leg-u1.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/boson-u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include <boost/program_options.hpp>

//
// Hamiltonian in original gauge (equation 1 of draft paper):
//
// H = -J \sum_r \sum_{l=0}^{2} a^\dagger_{l,r} a_{l,r+1} + h.c
//   + -J_\perp \sum_r \sum_{l=0}^{1} e^{-ir\phi} a^\dagger_{l,r} a_{l+1,r} + h.c.
// + U/2 \sum_{l,r} n_{l,r} (n_{l,r}-1)
//
// Additional term to Hamiltonian: cylindrical boundaries with flux theta.
// H_theta = e^{i\theta} a^\dagger_{2} a_{0} - e^{i\theta} a_{2} a^\dagger_{0}
//
// Change gauge to translationally invariant along the legs:
//
// b^\dagger_{0,r} = exp[-ir\phi] a^\dagger_{0,r}
// b^\dagger_{1,r} = a^\dagger_{1,r}
// b^\dagger_{2,r} = exp[ir\phi] a^\dagger_{0,r}
//
// New Hamiltonian:
// -J \sum_{r} ( exp[-i\phi] b^\dagger_{0,r+1} b_{0,r}
//                  + b^\dagger_{1,r+1} b_{1,r}
//                  + exp[i\phi] b^\dagger_{2,r+1} b_{2,r} + h.c.)
// -J_\perp \sum{r} b^\dagger_{0,r} b_{1,r} + b^\dagger_{1,r} b_{2,r} + H.c.
// + U/2 \sum_{l,r} n_{l,r} (n_{l,r}-1)
//
// For parameters J, K = J_\perp, U, alpha, the Hamiltonian is:
// J*(H_J1 + cos(pi*alpha)*(H_J2 + H_J0) + sin(pi*alpha)*(H_Jc2 - H_Jc0)) + K*H_K + U*H_U
//
// Current operator is derivative with respect to flux, so
//
// 
// (following equations are for 2-leg ladder!)
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
	 ("H_J2",  "nearest neighbor leg 2 hopping")
	 ("H_J",   "nearest neighbor leg hopping = H_J0 + H_J1 + H_J2")
	 ("H_Jc0", "nearest neighbor leg 0 current")
	 ("H_Jc1", "nearest neighbor leg 1 current")
	 ("H_Jc2", "nearest neighbor leg 2 current")
	 ("H_K",   "nearest neighbor rung hopping")
	 ("H_Kc",  "nearest neighbor rung current")
	 ("H_Kp",  "periodic rung hopping")
	 ("H_Kcp", "periodic rung current")
	 ("H_U",   "on-site Coulomb repulsion N*(N-1)/2")
	 ;

      OpDescriptions.add_functions()
	 ("H", "Hamiltonian, parametized by K, alpha (flux), U, J")
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
      UnitCell Cell(repeat(Site, 3));
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2"),
	 j0(Cell, "j0"), j1(Cell, "j1"), j2(Cell, "j2"), k0(Cell, "k0"), k1(Cell, "k1");

      j0 = std::complex<double>(0,1)*(BH(0)[0]*B(1)[0] - B(0)[0]*BH(1)[0]);
      j1 = std::complex<double>(0,1)*(BH(0)[1]*B(1)[1] - B(0)[1]*BH(1)[1]);
      j2 = std::complex<double>(0,1)*(BH(0)[2]*B(1)[2] - B(0)[2]*BH(1)[2]);

      k0 = std::complex<double>(0,1)*(BH(0)[0]*B(0)[1] - B(0)[0]*BH(0)[1]);
      k1 = std::complex<double>(0,1)*(BH(0)[1]*B(0)[2] - B(0)[1]*BH(0)[2]);

      InfiniteLattice Lattice(Cell);
      
      UnitCellMPO HJ0, HJ1, HJ2, HJc0, HJc1, HJc2, HK, HKc, HKp, HKcp, HU;
      HJ0 = BH(0)[0]*B(1)[0] + B(0)[0]*BH(1)[0];
      HJ1 = BH(0)[1]*B(1)[1] + B(0)[1]*BH(1)[1];
      HJ2 = BH(0)[2]*B(1)[2] + B(0)[2]*BH(1)[2];
      HJc0 = std::complex<double>(0,1)*(BH(0)[0]*B(1)[0] - B(0)[0]*BH(1)[0]);
      HJc1 = std::complex<double>(0,1)*(BH(0)[1]*B(1)[1] - B(0)[1]*BH(1)[1]);
      HJc2 = std::complex<double>(0,1)*(BH(0)[2]*B(1)[2] - B(0)[2]*BH(1)[2]);
      HU = 0.5*(N2(0)[0] + N2(0)[1] + N2(0)[2]);

      HK = BH(0)[0]*B(0)[1] + B(0)[0]*BH(0)[1] + BH(0)[1]*B(0)[2] + B(0)[1]*BH(0)[2];
      HKc = std::complex<double>(0,1)*(BH(0)[0]*B(0)[1] - B(0)[0]*BH(0)[1] + BH(0)[1]*B(0)[2] - B(0)[1]*BH(0)[2]);

      HKp = BH(0)[2]*B(0)[0] + B(0)[2]*BH(0)[0];
      HKcp = std::complex<double>(0,1)*(BH(0)[2]*B(0)[0] - B(0)[2]*BH(0)[0]);
      
      Lattice["H_J0"] = sum_unit(HJ0);
      Lattice["H_J1"] = sum_unit(HJ1);
      Lattice["H_J2"] = sum_unit(HJ2);
      Lattice["H_J"] = sum_unit(HJ0 + HJ1 + HJ2);
      Lattice["H_Jc0"] = sum_unit(j0);
      Lattice["H_Jc1"] = sum_unit(j1);
      Lattice["H_Jc2"] = sum_unit(j2);
      Lattice["H_K"] = sum_unit(HK);
      Lattice["H_Kc"] = sum_unit(HKc);
      Lattice["H_Kp"] = sum_unit(HKp);
      Lattice["H_Kcp"] = sum_unit(HKcp);
      Lattice["H_U"] = sum_unit(HU);

      // For parameters J, K = J_\perp, U, alpha, the Hamiltonian is:
      // J*(H_J1 + cos(pi*alpha)*(H_J2 + H_J0) + sin(pi*alpha)*(H_Jc2 - H_Jc0)) + K*H_K + U*H_U

      Lattice.func("H")("K", arg("alpha")=0.0, arg("U")=0.0, arg("J")=1.0)
	 = "J*(H_J1 + cos(pi*alpha)*(H_J2 + H_J0) + sin(pi*alpha)*(H_Jc2 - H_Jc0)) + K*H_K + U*H_U";

      Lattice.set_description("Bosonic 3-leg ladder with flux");
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
