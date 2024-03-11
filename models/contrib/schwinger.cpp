// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/schwinger.cpp
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
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
#include "models/spinlessfermion-u1.h"
#include "models/contrib/spinlessantifermion-u1-schwinger.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.description("Fermionic Schwinger model");
      OpDescriptions.author("Nicholas Godfrey", "n.godfrey@uq.net.au");
      OpDescriptions.add_operators()
         ("H_c"   , "Fermion hopping")
         ("H_m"   , "Fermion mass")
         ("H_ms"  , "Fermion mass, using 2*sum(phi^\\dagger phi)")
		   ("H_E"   , "Electric field (zero background)")
         ("H_Esym", "Electric field (symmetric form)")
         ;

      OpDescriptions.add_cell_operators()
         ("Nf" , "Number of fermions in a unit cell (+1 for fermion, -1 for antifermion)")
         ("Np" , "Number of particles in a unit cell (+1 for fermion, +1 for antifermion)")
         ;

      OpDescriptions.add_functions()
         ("H_l", "Background field l")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         std::cerr << "Write something about the prefactors in terms of g, a, m\n";
         return 1;
      }

      LatticeSite SiteEven = SpinlessFermionU1();
      LatticeSite SiteOdd = SpinlessAntifermionU1();

      UnitCell Cell(join(SiteEven, SiteOdd));

      InfiniteLattice Lattice(&Cell);
      UnitCellOperator I(Cell, "I"), CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N"), Nf(Cell, "Nf"), Np(Cell, "Np"), P(Cell, "P");

      // Short-cut operators on the unit cell
      Nf = Nf(0)[0] + Nf(0)[1];
      Np = Np(0)[0] + Np(0)[1];

      Lattice["H_c"] = std::complex<double>(0,-1) * sum_unit
         (
              dot(CH(0)[0], C(0)[1])
            + dot(C(0)[0], CH(0)[1])
            + dot(CH(0)[1], C(1)[0])
            + dot(C(0)[1], CH(1)[0])
         );

      Lattice["H_m"] = 2*sum_unit(Np(0));
      Lattice["H_ms"] = 2*sum_unit(Np(0) - I(0));

      Lattice["H_E"] = sum_partial(2 * pow(sum_unit(Nf(0)),2),I)
         + 2 * sum_partial(sum_unit(Nf(0)),N(0)[0])
         + sum_unit(pow(N(0)[0],2));  // the pow here is redundant since N(0)[0] is a projector, but it is technically there
         //+ pow(sum_unit(Nf(0)),2);

      // Symmeterized version of H_E
      Lattice["H_Esym"] = sum_string(-Nf(0), I(0)[0], 2*I(0)[0], I(0)[0], Nf(0))
         + sum_string(-Nf(0), I(0)[0], Nf(0)[1])
         + sum_string(-Nf(0) - N(0)[0], I(0)[0], Nf(0))
         + sum_unit(-Nf(0)[0]*Nf(0)[1]);

      Lattice.func("H_l")("l") = "l*2*(2 * sum_partial(sum_unit(Nf(0)),I) + sum_unit(2*Np(0)[0] - Np(0)[1]))";

      // Information about the lattice
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice
      pheap::ExportObject(FileName, Lattice);
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
