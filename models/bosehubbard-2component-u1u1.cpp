// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/bosehubbard-2component-u1u1.cpp
//
// Copyright (C) 2015-2023 Ian McCulloch <ian@qusim.net>
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
#include "models/boson-2component-u1u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include <boost/program_options.hpp>

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
      OpDescriptions.set_description("Bosonic 2-leg ladder with U(1)xU(1) symmetry");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_JA"   , "nearest-neighbor hopping for the A-species")
         ("H_JB"   , "nearest-neighbor hopping for the B-species")
         ("H_UA"   , "Coulomb repulsion for species A")
         ("H_UB"   , "Coulomb repulsion for species B")
         ("H_U"    , "Intra-species coulomb repulsion, H_UA + H_UB")
         ("H_UAB"  , "inter-species Coulomb repulsion")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << 'n';
         return 1;
      }

      LatticeSite Site = Boson2ComponentU1U1(MaxN);
      UnitCell Cell = Site;
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator BH_A(Cell, "BH_A"), B_A(Cell, "B_A"), N_A(Cell, "N_A"), N2_A(Cell, "N2_A"),
         BH_B(Cell, "BH_B"), B_B(Cell, "B_B"), N_B(Cell, "N_B"), N2_B(Cell, "N2_B"), I(Cell, "I");

      Lattice["H_JA"]   = -sum_unit(BH_A(0)*B_A(1) + B_A(0)*BH_A(1));
      Lattice["H_JB"]   = -sum_unit(BH_B(0)*B_B(1) + B_B(0)*BH_B(1));
      Lattice["H_UA"]   = 0.5*sum_unit(N_A(0)*(N_A(0)-I(0)));
      Lattice["H_UB"]   = 0.5*sum_unit(N_B(0)*(N_B(0)-I(0)));
      Lattice["H_U"]    = Lattice["H_UA"] + Lattice["H_UB"];
      Lattice["H_UAB"]  = sum_unit(N_A(0)*N_B(0));

      // Add projectors onto N+1 and N-1 bosons

      UnitCellOperator NP(Cell, "NP"), NM(Cell, "NM");
      LatticeSite NPSite = Boson2ComponentU1U1(MaxN+1);
      LatticeSite NMSite = Boson2ComponentU1U1(MaxN-1);
      SimpleOperator NPOp(NPSite.Basis1(), Site.Basis1());
      SimpleOperator NMOp(NMSite.Basis1(), Site.Basis1());
      for (int na = 0; na <= MaxN; ++na)
      {
         for (int nb = 0; nb <= MaxN; ++nb)
         {
            int n = na*(MaxN+1)+nb;
            int np = na*(MaxN+2)+nb;
            NPOp(np, n) = 1.0;
            if (na < MaxN && nb < MaxN)
            {
               int nm = na*MaxN+nb;
               NMOp(nm, n) = 1.0;
            }
         }
      }
      NP(0) = Cell.map_local_operator(NPOp, LatticeCommute::Bosonic, "Projector onto N+1 bosons", 0, 0);
      NM(0) = Cell.map_local_operator(NMOp, LatticeCommute::Bosonic, "Projector onto N-1 bosons", 0, 0);

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
