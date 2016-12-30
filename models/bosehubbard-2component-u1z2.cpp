// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/bosehubbard-2component-u1z2.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "models/boson-2component-u1z2.h"
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
      OpDescriptions.set_description("Bosonic 2-leg ladder with flux");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J"   , "nearest-neighbor hopping")
         ("H_K"   , "tunnelling  between components")
         ("H_U"   , "intra-species Coulomb repulsion")
         ("H_U12" , "inter-species Coulomb repulsion")
         ("D"     , "difference in occupation number between components\n")
         ("D2"    , "squared difference in occuptation number between components\n")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << 'n';
         return 1;
      }

      LatticeSite Site = Boson2ComponentU1Z2(MaxN);
      UnitCell Cell = Site;
      InfiniteLattice Lattice(&Cell);

      UnitCellOperator BH_A(Cell, "BH_A"), B_A(Cell, "B_A"), N_A(Cell, "N_A"), N2_A(Cell, "N2_A"),
         BH_S(Cell, "BH_S"), B_S(Cell, "B_S"), N_S(Cell, "N_S"), N2_S(Cell, "N2_S");


      UnitCellMPO HJ = -(BH_A(0)*B_A(1) + B_A(0)*BH_A(1) + BH_S(0)*B_S(1) + B_S(0)*BH_S(1));
      UnitCellMPO HK = -(N_S(0) - N_A(0));

      UnitCellMPO PairHopping = pow(BH_S(0)*B_A(0),2) + pow(BH_A(0)*B_S(0),2);

      UnitCellMPO HU = N_S(0)*N_A(0) + 0.25 * (N2_S(0) + N2_A(0) + PairHopping);
      UnitCellMPO HU12 = 0.25 * (N2_S(0) + N2_A(0) - PairHopping);

      UnitCellMPO D = BH_A(0)*B_S(0) + BH_S(0)*B_A(0);  // the order parameter, is antisymmetric in this basis

      Lattice["H_J"]   = sum_unit(HJ);
      Lattice["H_K"]   = sum_unit(HK);
      Lattice["H_U"]   = sum_unit(HU);
      Lattice["H_U12"] = sum_unit(HU12);
      Lattice["D"]     = sum_unit(D);
      Lattice["D2"]    = sum_unit(D*D);

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
