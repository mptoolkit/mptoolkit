// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/u1-qlm-u1.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2020-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// U(1) 1+1D U(1) quantum link model

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spinlessfermion-u1.h"
#include "models/spinlessantifermion-u1.h"
#include "models/spin.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      half_int Spin = 0.5;
      bool Tau = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), FormatDefault("magnitude of the link spins", Spin).c_str())
         ("tau", prog_opt::bool_switch(&Tau), "use alternative coeffients for the ladder operators based on the boson ladder operators")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) 1+1D U(1) quantum link model");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_cell_operators()
         ("G0"   , "Gauss's law operator for the fermion site 0")
         ("G2"   , "Gauss's law operator for the antifermion site 2")
         ;
      OpDescriptions.add_operators()
         ("H_t"  , "nearest-neighbor hopping")
         ("H_m"  , "fermion mass")
         ("H_chi", "background field")
         ("H_g"  , "gauge coupling")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         std::cerr << "Operators:" << std::endl << OpDescriptions;
         return 1;
      }

      LatticeSite FSite = SpinlessFermionU1();
      LatticeSite AFSite = SpinlessAntifermionU1();
      LatticeSite SSite1 = SpinSite(Spin);
      LatticeSite SSite2 = SpinSite(Spin);

      if (Tau)
      {
         std::map<half_int, std::string> SpinBasis;
         for (half_int s = -Spin; s <= Spin; ++s)
            SpinBasis[s] = boost::lexical_cast<std::string>(s);

         for (half_int s = -Spin; s < Spin; ++s)
         {
            int n = (Spin + s).twice();
            double Coeff = std::sqrt((n + 1) * (n + 2));
            SSite1["Sp"](SpinBasis[-s], SpinBasis[-s-1]) = Coeff;
            SSite2["Sp"](SpinBasis[s+1], SpinBasis[s]) = Coeff;
         }

         SSite1["Sm"] = adjoint(SSite1["Sp"]);
         SSite2["Sm"] = adjoint(SSite2["Sp"]);
      }

      UnitCell FCell(FSite.GetSymmetryList(), FSite, SSite1);
      UnitCell AFCell(AFSite.GetSymmetryList(), AFSite, SSite2);
      UnitCell Cell(join(FCell, AFCell));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), N(Cell, "N"),
                       Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");

      UnitCellMPO t, m, chi, g;

      t += Sp(0)[1] * dot(CH(0)[0], C(0)[2]) + Sm(0)[1] * dot(CH(0)[2], C(0)[0]);
      t += Sp(0)[3] * dot(CH(0)[2], C(1)[0]) + Sm(0)[3] * dot(CH(1)[0], C(0)[2]);

      m += N(0)[0] - N(0)[2];

      chi += Sz(0)[1] + Sz(0)[3];
      g += Sz(0)[1]*Sz(0)[1] + Sz(0)[3]*Sz(0)[3];

      Lattice["H_t"] = sum_unit(t);
      Lattice["H_m"] = sum_unit(m);
      Lattice["H_chi"] = sum_unit(chi);
      Lattice["H_g"] = sum_unit(g);

      // Gauss's law operators.
      UnitCellOperator G0(Cell, "G0"), G2(Cell, "G2");

      G0 = N(0)[0] - Sz(0)[1] + Sz(-1)[3];
      G2 = N(0)[2] - Sz(0)[3] + Sz(0)[1];

      // Information about the lattice
      Lattice.set_command_line(argc, argv);
      Lattice.set_operator_descriptions(OpDescriptions);

      // save the lattice to disk
      pheap::ExportObject(FileName, Lattice);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl; 
      return 1;
   }
}
