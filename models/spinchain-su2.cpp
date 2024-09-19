// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/spinchain-su2.cpp
//
// Copyright (C) 2004-2019 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2015-2016 Seyed Saadatmand <s.saadatmand@uq.edu.au>
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
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

// spin-1/2 chain exact energy per site is 1/4 - ln(2)

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      // Parameters of the lattice (with defaults, if applicable)
      half_int Spin = 0.5;
      std::string FileName;

      // Options
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

      // Descriptions of each operator
      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("SU(2) spin chain");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1"  , "nearest neighbor spin exchange")
         ("H_J2"  , "next-nearest neighbor spin exchange")
         ("H_J3"  , "next-next-nearest neighbor spin exchange")
         ("H_D"   , "-3*P_D, projector onto singlet dimer")
         ("H_T"   , "-6*P_T, projector onto singlet trimer")
         ("H_B1"  , "nearest neighbor biquadratic spin exchange (S.S)^2")
         ("H_B2"  , "next-nearest neighbor biquadratic spin exchange (S.S)^2")
         ("H_B3"  , "next-next-nearest neighbor biquadratic spin exchange (S.S)^2")
         ("H_Q1"  , "nearest neighbor quadrupole exchange (Q.Q)")
         ("H_Q2"  , "next-nearest neighbor quadrupole exchange (Q.Q)")
         ("H_Q3"  , "next-next-nearest neighbor quadrupole exchange (Q.Q)")
         ("H_AKLT", "AKLT Hamiltonian H_J1 + (1/3)*H_J2", "spin 1", [&Spin]()->bool {return Spin==1;})
         ;

      // Descriptions for the operators
      OpDescriptions.add_functions()
         ("H_exp", "Exponential decay spin exchange parameterized by lambda as exp(-lambda*r)")
	 ("H_BQ"  , "Bilinear-biquadratic model, parameterized by theta", "spin 1", [&Spin]()->bool {return Spin==1;})
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite Site = SpinSU2(Spin);

      // The UnitCell consists of a single site
      UnitCell Cell(Site);

      // Make an infinite lattice of our unit cell
      InfiniteLattice Lattice(&Cell);

      // A short-cut to refer to an operator defined within our unit cell
      UnitCellOperator S(Cell, "S"), Q(Cell, "Q"), I(Cell, "I");

      // Define operators that have support over the infinite lattice
      Lattice["H_J1"] = sum_unit(inner(S(0), S(1)));
      Lattice["H_J2"] = sum_unit(inner(S(0), S(2)));
      Lattice["H_J3"] = sum_unit(inner(S(0), S(3)));

      Lattice["H_B1"] = sum_unit(pow(inner(S(0), S(1)), 2));
      Lattice["H_B2"] = sum_unit(pow(inner(S(0), S(2)), 2));
      Lattice["H_B3"] = sum_unit(pow(inner(S(0), S(3)), 2));

      Lattice["H_Q1"] = sum_unit(inner(Q(0), Q(1)));
      Lattice["H_Q2"] = sum_unit(inner(Q(0), Q(2)));
      Lattice["H_Q3"] = sum_unit(inner(Q(0), Q(3)));

      Lattice["H_D"] = -sum_unit(pow(inner(S(0), S(1)), 2) - I(0));
      UnitCellMPO S3 = inner(S(0)+S(1)+S(2), S(0)+S(1)+S(2));
      Lattice["H_T"] = (1.0/24.0) * sum_unit((S3 - 2*I(0)) * (S3 - 6*I(0)) * (S3 - 12*I(0)));

      if (Spin == 1)
      {
         Lattice["H_AKLT"] = Lattice["H_J1"] + (1.0/3.0)*Lattice["H_B1"];
	 Lattice.func("H_BQ")("theta") = "cos(theta)*H_J1 + sin(theta)*H_B1";
      }

      Lattice.func("H_exp")(arg("lambda") = 0.5)
                  = "exp(-lambda)*sum_string_inner( S(0), exp(-lambda)*I(0), S(0) )";

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
