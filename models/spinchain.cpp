// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinchain.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
// Ising exact energy 4/pi per site.
// finite size OBC: 1 / sin(pi / (4L+w))

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
      OpDescriptions.set_description("Spin chain");
      OpDescriptions.add_operators()
	 ("H_xx", "nearest neighbor spin coupling Sx Sx")
	 ("H_yy", "nearest neighbor spin exchange Sy Sy")
	 ("H_zz", "nearest neighbor spin exchange Sz Sz")
	 ("H_x" , "magnetic field in the x direction")
	 ("H_y" , "magnetic field in the y direction")
	 ("H_z" , "magnetic field in the z direction")
	 ("H_J1z", "same as H_zz")
	 ("H_J1t", "transverse spin exchange, H_xx + H_yy")
	 ("H_J1" , "nearest neighbor spin exchange = H_J1z + H_J1t")
	 ("H_B1" , "nearest neighbor biquadratic spin exchange (S.S)^2")
	 ("H_mu" , "single-ion anistotropy, H_mu = sum_i Sz(i)^2")
	 ("X"    , "pi rotation about the X axis, equivalent to prod_unit(exp(i*pi*Sx(0)))")
	 ("Y"    , "pi rotation about the Y axis, equivalent to prod_unit(exp(i*pi*Sy(0)))")
	 ("Z"    , "pi rotation about the Z axis, equivalent to prod_unit(exp(i*pi*Sz(0)))")
	 ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << OpDescriptions << '\n';
	 std::cerr << "only for spin-1: H_AKLT  - AKLT Hamiltonian H_J1 + (1/3)*H_B1\n";
         return 1;
      }

      LatticeSite Site = SpinSite(Spin);
      UnitCell Cell(Site);
      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"); // identity operator
      UnitCellOperator Test(Cell, "Test");

      Test(0) = 0.0 * I(0);

      for (int i = 0; i < 50; ++i)
      {
	 Test(0) = Test(0) + 0.5 * Sx(i)*Sx(i+1);
	 Test(0) = Test(0) + (1.0/3.0) * Sy(i)*Sy(i+1);
	 Test(0) = Test(0) + 0.2 * Sz(i)*Sz(i+1);
      }

      InfiniteLattice Lattice("Spin chain", Cell);

      UnitCellMPO SpinExchange = Sx(0)*Sx(1) + Sy(0)*Sy(1) + Sz(0)*Sz(1);

      Lattice["H_xx"] = sum_unit(Sx(0)*Sx(1));
      Lattice["H_yy"] = sum_unit(Sy(0)*Sy(1));
      Lattice["H_zz"] = sum_unit(Sz(0)*Sz(1));

      Lattice["H_x"] = sum_unit(Sx(0));
      Lattice["H_y"] = sum_unit(Sy(0));
      Lattice["H_z"] = sum_unit(Sz(0));

      Lattice["H_J1z"] = Lattice["H_zz"];
      Lattice["H_J1t"] = Lattice["H_xx"] + Lattice["H_yy"];
      Lattice["H_J1"] = sum_unit(SpinExchange);
      Lattice["H_B1"] = sum_unit(SpinExchange*SpinExchange);

      Lattice["H_mu"] = sum_unit(Sz(0)*Sz(0));

      if (Spin == 1)
      {
	 Lattice["H_AKLT"] = Lattice["H_J1"] + (1.0/3.0)*Lattice["H_B1"];
	 Lattice["H_AKLT"].set_description("AKLT Hamiltonian H_J1 + (1/3)*H_B1");
      }

      Lattice["X"] = prod_unit(exp(std::complex<double>(0.0, math_const::pi)*Sx(0)));
      Lattice["Y"] = prod_unit(exp(std::complex<double>(0.0, math_const::pi)*Sy(0)));
      Lattice["Z"] = prod_unit(exp(std::complex<double>(0.0, math_const::pi)*Sz(0)));

      /* Lattice.func("HaldShast")(arg("lambda") = 0.5)
                  = "sum_kink( (1/lambda)*I(0), Sz(0) ) * sum_kink( lambda*I(0), Sz(0) ) + sum_kink( (1/lambda)*I(0), Sy(0) ) * sum_kink( lambda*I(0), Sy(0) ) + sum_kink( (1/lambda)*I(0), Sx(0) ) * sum_kink( lambda*I(0), Sx(0) )"; */

      /* Lattice.func("Test")(arg("J") = 0.0)
         = "J*sum_unit( I(0)*Sx(0)*Sx(1) )"; */

      // Information about the lattice
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
