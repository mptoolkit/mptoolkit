// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/bosehubbardcylinder-u1.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "models/boson-u1.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      int MaxN = DefaultMaxN;
      int x = 0;
      int y = 4;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&MaxN),
          FormatDefault("Maximum number of bosons per site", MaxN).c_str())
         (",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str())
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector", y).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) Bose-Hubbard cylinder square lattice");
      OpDescriptions.author("J Osborne", "j.osborne@uqconnect.edu.au");
      OpDescriptions.add_operators()
         ("H_Jxp"  , "nearest-neighbor hopping in +x-direction")
         ("H_Jxm"  , "nearest-neighbor hopping in -x-direction")
         ("H_Jyp"  , "nearest-neighbor hopping in +y-direction")
         ("H_Jym"  , "nearest-neighbor hopping in -y-direction")
         ("H_Jx"   , "nearest-neighbor hopping in x-direction")
         ("H_Jy"   , "nearest-neighbor hopping in y-direction")
         ("H_J"    , "nearest-neighbor hopping")
         ("H_U"    , "on-site Coulomb interaction N*(N-1)/2")
         ("H_Vx"   , "nearest-neighbour Coulomb repulsion in x-direction")
         ("H_Vy"   , "nearest-neighbour Coulomb repulsion in x-direction")
         ("H_V"    , "nearest-neighbour Coulomb repulsion")
         ("H_delta", "QLM link potential", "x = 0, y even",
         [&x, &y]()->bool{return x == 0 && y%2 == 0;})
         ("H_eta"  , "QLM eta potential", "x = 0, y even",
         [&x, &y]()->bool{return x == 0 && y%2 == 0;})
         ;
      OpDescriptions.add_functions()
         ("H_Jx"   , "nearest-neighbor complex hopping in x-direction")
         ("H_Jy"   , "nearest-neighbor complex hopping in x-direction")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      int CellSize = x == 0 ? y : x;

      LatticeSite Site = BosonU1(MaxN);
      UnitCell Cell(repeat(Site, CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2");

      // Define hopping terms and near-neighbour interactions.
      UnitCellMPO Jxp, Jxm, Jyp, Jym, Vx, Vy;

      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            Jxp += BH(0)[i]*B(1)[i];
            Jxm += B(0)[i]*BH(1)[i];
            Jyp += BH(0)[i]*B(0)[(i+1)%y];
            Jym += B(0)[i]*BH(0)[(i+1)%y];
            Vx += N(0)[i]*N(1)[i];
            Vy += N(0)[i]*N(0)[(i+1)%y];
         }
      }
      else
      {
         for (int i = 0; i < x-1; ++i)
         {
            Jxp += BH(0)[i]*B(0)[i+1];
            Jxm += B(0)[i]*BH(0)[i+1];
            Vx += N(0)[i]*N(0)[i+1];
         }
         Jxp += BH(0)[x-1]*B(y+1)[0];
         Jxm += B(0)[x-1]*BH(y+1)[0];
         Vx += N(0)[x-1]*N(y+1)[0];

         for (int i = 0; i < x; ++i)
         {
            Jyp += BH(0)[i]*B(1)[i];
            Jym += B(0)[i]*BH(1)[i];
            Vy += N(0)[i]*N(1)[i];
         }
      }

      Lattice["H_Jxp"] = sum_unit(Jxp);
      Lattice["H_Jxm"] = sum_unit(Jxm);
      Lattice["H_Jx"] = sum_unit(Jxp+Jxm);
      Lattice.func("H_Jx")(arg("theta")) = "exp(-i*theta)*H_Jxp + exp(i*theta)*H_Jxm";
      Lattice["H_Jyp"] = sum_unit(Jyp);
      Lattice["H_Jym"] = sum_unit(Jym);
      Lattice["H_Jy"] = sum_unit(Jyp+Jym);
      Lattice.func("H_Jy")(arg("theta")) = "exp(-i*theta)*H_Jyp + exp(i*theta)*H_Jym";
      Lattice["H_J"] = sum_unit(Jxp+Jxm+Jyp+Jym);

      Lattice["H_Vx"] = sum_unit(Vx);
      Lattice["H_Vy"] = sum_unit(Vy);
      Lattice["H_V"] = sum_unit(Vx+Vy);

      // Define on-site interactions.
      UnitCellMPO U;

      for (int i = 0; i < CellSize; ++i)
         U += 0.5*N2(0)[i];

      Lattice["H_U"] = sum_unit(U);

      // Define QLM potentials.
      if (x == 0 && y%2 == 0)
      {
         UnitCellMPO H_delta, H_eta;

         for (int i = 0; i < CellSize/2; ++i)
         {
            H_delta += N(0)[2*i+1] + N(1)[2*i];
            H_eta += N(1)[2*i+1];
         }

         Lattice["H_delta"] = sum_unit(H_delta, 2*CellSize);
         Lattice["H_eta"] = sum_unit(H_eta, 2*CellSize);
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
