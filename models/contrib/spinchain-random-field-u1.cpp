// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spinchain-random-field-u1.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "models/spin-u1.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/randutil.h"
#include <sstream>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      std::string FileName;
      int Length = 0;
      uint32_t Seed = 0;
      std::string FieldStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
	 ("unitcell,u", prog_opt::value(&Length), "unit cell size")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
	 ("seed,s", prog_opt::value(&Seed), "random seed [if not specified, this is initialized randomly]")
         ("field", prog_opt::value(&FieldStr), "user-specified random field configuration into H_Field")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style).
                      //                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) spin chain with random field");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1z"  , "nearest neighbor spin coupling Sz Sz")
         ("H_J1t"  , "nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
         ("H_J1"   , "nearest neighbor spin exchange = H_J1z + H_J1t")
         ("H_Field", "user-specified field configuration", "option --field", [&FieldStr]()->bool {return !FieldStr.empty();})
         ;
      OpDescriptions.add_functions()
         ("Field" , "random bonds uniform on [-1,1], indexed by an integer disorder realization")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      if (!vm.count("seed"))
      {
         // we only need 1 random number, so might as well use crypto_rand(),
         // since this consumes fewer resources than seeding the RNG
         // seed randomly
         Seed = randutil::crypto_rand();
      }

      std::vector<double> Field;
      std::istringstream Str(FieldStr);
      double x;
      while (Str >> x)
      {
         Field.push_back(x);
      }

      if (Length == 0)
         Length = Field.size();

      if (Length == 0)
      {
         std::cerr << basename(argv[0]) << ": fatal: unit cell length must be specified.\n";
         return 1;
      }

      std::string FieldFunc = "sum_unit((2*rand(SEED,n,0)-1)*Sz(0)[0]";
      for (int i = 1; i < Length; ++i)
      {
         FieldFunc += "+(2*rand(SEED,n," + std::to_string(i) + ")-1)*Sz(0)[" + std::to_string(i) + "]";
      }
      FieldFunc += ')';

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(repeat(Site, Length));
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      InfiniteLattice Lattice(&Cell);

      if (!Field.empty())
      {
         if (Field.size() != Length)
         {
            std::cerr <<  basename(argv[0]) << ": fatal: --fields option must supply name number of terms as the unit cell size.\n";
            return 1;
         }
         UnitCellMPO F;
         for (int i = 0; i < Field.size(); ++i)
         {
            F += Field[i] * Sz(0)[i];
         }
         Lattice["H_Field"] = sum_unit(F);
      }

      for (int i = 0; i < Length-1; ++i)
      {
	 Lattice["H_J1z"] += sum_unit(Sz(0)[i]*Sz(0)[i+1]);
	 Lattice["H_J1t"] += 0.5 * sum_unit(Sp(0)[i]*Sm(0)[i+1] + Sm(0)[i]*Sp(0)[i+1]);
      }
      Lattice["H_J1"] = Lattice["H_J1z"] + Lattice["H_J1t"];

      Lattice.func("Field")(arg("n")) = FieldFunc;

      Lattice.arg("SEED") = std::complex<double>(Seed);

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
