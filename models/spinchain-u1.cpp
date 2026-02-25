// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spinchain-u1.cpp
//
// Copyright (C) 2004-2023 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

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
      OpDescriptions.set_description("U(1) Spin chain");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1z"  , "nearest neighbor spin coupling Sz Sz")
         ("H_J1t"  , "nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
         ("H_J1"   , "nearest neighbor spin exchange = H_J1z + H_J1t")
         ("H_J2z"  , "next-nearest neighbor spin coupling Sz Sz")
         ("H_J2t"  , "next-nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
         ("H_J2"   , "next-nearest neighbor spin exchange = H_J1z + H_J1t")
         ("H_B1"   , "nearest neighbor biquadratic spin exchange (S.S)^2")
         ("H_B2"   , "next-nearest neighbor biquadratic spin exchange (S.S)^2")
         ("H_B1xy" , "nearest neighbor biquadratic XY spin exchange (Sx.Sx + Sy.Sy)^2")
         ("H_mu"   , "single-ion anistotropy, H_mu = sum_i Sz(i)^2")
         ("H_PD"   , "spin-1 projector onto singlet dimers (S(0)+S(1))^2 - 2) * (S(0)+S(1))^2 - 6)")
         ("H_PT"   , "spin-1 projector onto singlet trimers (S(0)+S(1)+S(2))^2 - 2)*(S(0)+S(1)+S(2))^2 - 6)*(S(0)+S(1)+S(2))^2 - 12)")
         ("H_dimer", "dimerized spin exchange, sum_i S(2*i).S(2*i+1) - S(2*i+1).S(2*i+2)")
         ("H_stag" , "staggered field (-1)^n Sz(n)")
         ("H_AKLT" , "AKLT Hamiltonian H_J1 + (1/3)*H_B1", "spin 1", [&Spin]()->bool {return Spin==1;})
         ;

      OpDescriptions.add_functions()
	 ("H_BQ"  , "Bilinear-biquadratic model, parameterized by theta", "spin 1",
	  [&Spin]()->bool {return Spin==1;})
         ("H_murray", "Biquadratic model with anisotropy, parameterized by xy and z", "spin 1",
	  [&Spin]()->bool {return Spin==1;})
	 ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(Site);
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"), SS(Cell, "SS"), S_2(Cell, "S_2"), S_3(Cell, "S_3");
      InfiniteLattice Lattice(&Cell);

      // Nearest neighbor S(0).S(1)
      SS = Sz(0)*Sz(1)+ 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1));
      SS.set_description("Shortcut for S(0).S(1)");
      // (S(0)+S(1))^2
      S_2 = 4*I(0) + 2*SS(0);
      S_2.set_description("Shortcut for (S(0)+S(1))^2 (spin-1 only)");
      // (S(0)+S(1)+S(2))^2
      S_3 = 6*I(0) + 2*(SS(0) + SS(1) + Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));
      S_3.set_description("Shortcut for (S(0)+S(1)+S(2))^2 (spin-1 only)");

      Lattice["H_J1z"] = sum_unit(Sz(0)*Sz(1));
      Lattice["H_J1t"] = 0.5 * sum_unit(Sp(0)*Sm(1) + Sm(0)*Sp(1));
      Lattice["H_J1"]  = Lattice["H_J1z"] + Lattice["H_J1t"];

      Lattice["H_J2z"] = sum_unit(Sz(0)*Sz(2));
      Lattice["H_J2t"] = sum_unit(0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));
      Lattice["H_J2"] = sum_unit(Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)));

      Lattice["H_B1"] = sum_unit(pow(Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)), 2));
      Lattice["H_B2"] = sum_unit(pow(Sz(0)*Sz(2) + 0.5*(Sp(0)*Sm(2) + Sm(0)*Sp(2)), 2));

      Lattice["H_B1xy"] = sum_unit(pow(0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)), 2));

      Lattice["H_mu"] = sum_unit(Sz(0)*Sz(0));

      Lattice["H_dimer"] = sum_unit(Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1))
				    - (Sz(1)*Sz(2) + 0.5*(Sp(1)*Sm(2) + Sm(1)*Sp(2))), 2);

      Lattice["H_stag"] = sum_unit(Sz(0) - Sz(1), 2);

      Lattice["H_PD"] = (1.0 / 12.0) * sum_unit( (S_2(0) - 2*I(0)) * (S_2(0) - 6*I(0)) );
      Lattice["H_PT"] = (-1.0 / 144.0) * sum_unit( (S_3(0) - 2*I(0)) * (S_3(0) - 6*I(0)) * (S_3(0) - 12*I(0)) );


      if (Spin == 1)
      {
         Lattice["H_AKLT"] = Lattice["H_J1"] + (1.0/3.0)*Lattice["H_B1"];
	 Lattice.func("H_BQ")("theta") = "cos(theta)*H_J1 + sin(theta)*H_B1";
         // note the parser doesn't handle precedence of negation and power properly,
         // -3^2 == (-3)^2 rather than -(3^2).  So we need to include brackets.
         Lattice.func("H_murray")("xy", "z") = "sum_unit(-((xy*0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1)) "
            " + z*Sz(0)*Sz(1))^2))";
      }

      Lattice.func("H_J1t_twist")("theta") =
         "0.5 * sum_unit(exp(i*theta)*Sp(0)*Sm(1) + exp(-i*theta)*Sm(0)*Sp(1))";
      Lattice.func("H_J1t_twist").set_description("Twisted spin exchange");

      Lattice.func("H_J1_twist")("theta") = "H_J1z + H_J1t_twist{theta}";
      Lattice.func("H_J1_twist").set_description("Twisted J1");

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
