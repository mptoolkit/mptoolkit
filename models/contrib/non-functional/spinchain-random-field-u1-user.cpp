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
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>


namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      std::string FileName;
      std::vector<double> Field;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("op", prog_opt::value(&OpStr), "op")
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("op", 1);
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) spin chain with random field");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J1z" , "nearest neighbor spin coupling Sz Sz")
         ("H_J1t" , "nearest neighbor spin exchange (1/2)(Sp Sm + Sm Sp)")
         ("H_J1"  , "nearest neighbor spin exchange = H_J1z + H_J1t")
	 ("H_h"   , "user-specified field configuration")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions;
         return 1;
      }

      std::vector<double> Fields = {0.781016, 1.48551, 3.44303, 5.50825, 1.64421, 5.72713, 0.718131, 5.55603,
				    -1.22152, 1.97702, 2.33431, -1.84989, -0.998605, -3.23945, 6.26837, -7.09259 };
      TRACE(Fields.size());
      double h = 8;
      {
	 boost::mt19937 generator;
	 generator.seed(0);
	 boost::random::uniform_real_distribution<double> box(-h,+h);
	 {
	    std::cout << "# fields= {";
	    for(size_t i=0; i<Fields.size(); i++)
	    {
	       Fields[i]= box(generator);
	       std::cout << Fields[i] << " ";
	    }
	 }
	 std::cout << '}' << std::endl;
      }

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell(repeat(Site, Length));
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      InfiniteLattice Lattice(&Cell);

      for (int i = 0; i < Length-1; ++i)
      {
	 Lattice["H_J1z"] += sum_unit(Sz(0)[i]*Sz(0)[i+1]);
	 Lattice["H_J1t"] += 0.5 * sum_unit(Sp(0)[i]*Sm(0)[i+1] + Sm(0)[i]*Sp(0)[i+1]);
      }
      Lattice["H_J1"] = Lattice["H_J1z"] + Lattice["H_J1t"];

      for (int i = 0; i < Length; ++i)
      {
	 Lattice["H_h"] += sum_unit(Fields[i] * Sz(0)[i]);
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
