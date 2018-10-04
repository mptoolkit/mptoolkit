// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/spinchain-z2.cpp
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

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-z2.h"
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
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_xx", "nearest neighbor spin coupling Sx Sx")
         ("H_yy", "nearest neighbor spin exchange Sy Sy")
         ("H_zz", "nearest neighbor spin exchange Sz Sz")
         ("H_x" , "magnetic field in the x direction")
         ("H_J1z", "same as H_zz")
         ("H_J1t", "transverse spin exchange, H_xx + H_yy")
         ("H_J1" , "nearest neighbor spin exchange = H_J1z + H_J1t")
         ("H_B1" , "nearest neighbor biquadratic spin exchange (S.S)^2")
         ("H_mu" , "single-ion anistotropy, H_mu = sum_i Sz(i)^2")
         ("H_ITF", "Ising transverse-field model, equivalent to -4*H_zz + 2*H_x")
         ("H_AKLT" , "AKLT Hamiltonian H_J1 + (1/3)*H_B1", "spin 1", [&Spin]()->bool {return Spin==1;})
         ;

      OpDescriptions.add_functions()
	 ("H_BQ"  , "Bilinear-biquadratic model, parameterized by theta", "spin 1",
	  [&Spin]()->bool {return Spin==1;})
         ("H_murray", "Biquadratic model with anisotropy, parameterized by x, y and z", "spin 1")
          ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      LatticeSite Site = SpinZ2(Spin);
      UnitCell Cell(Site);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator Sx(Cell, "Sx"), Sy(Cell, "Sy"), Sz(Cell, "Sz"), Sp(Cell, "Sp"), Sm(Cell, "Sm");
      UnitCellOperator I(Cell, "I"); // identity operator

      UnitCellMPO SpinExchange = Sx(0)*Sx(1) + Sy(0)*Sy(1) + Sz(0)*Sz(1);

      Lattice["H_xx"] = sum_unit(Sx(0)*Sx(1));
      Lattice["H_yy"] = sum_unit(Sy(0)*Sy(1));
      Lattice["H_zz"] = sum_unit(Sz(0)*Sz(1));

      Lattice["H_x"] = sum_unit(Sx(0));

      Lattice["H_J1z"] = Lattice["H_zz"];
      Lattice["H_J1t"] = Lattice["H_xx"] + Lattice["H_yy"];
      Lattice["H_J1"] = sum_unit(SpinExchange);
      Lattice["H_B1"] = sum_unit(SpinExchange*SpinExchange);

      Lattice["H_mu"] = sum_unit(Sz(0)*Sz(0));

      Lattice["H_ITF"] = -4*Lattice["H_zz"] + 2*Lattice["H_x"];

      if (Spin == 1)
      {
         Lattice["H_AKLT"] = Lattice["H_J1"] + (1.0/3.0)*Lattice["H_B1"];
	 Lattice.func("H_BQ")("theta") = "cos(theta)*H_J1 + sin(theta)*H_B1";
         // note the parser doesn't handle precedence of negation and power properly,
         // -3^2 == (-3)^2 rather than -(3^2).  So we need to include brackets.
         Lattice.func("H_murray")("x", arg("y")=1, arg("z")=1) = "sum_unit(-((x*Sx(0)*Sx(1) + y*Sy(0)*Sy(1) + z*Sz(0)*Sz(1))^2))";
      }

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
