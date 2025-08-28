// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/bosehubbardcylinder-qlm-u1.cpp
//
// Copyright (C) 2025 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
      half_int QLMSpin = 0.5;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&MaxN),
          FormatDefault("maximum number of bosons per site", MaxN).c_str())
         (",x", prog_opt::value(&x), FormatDefault("x wrapping vector", x).c_str())
         (",y", prog_opt::value(&y), FormatDefault("y wrapping vector", y).c_str())
         ("qlm-spin", prog_opt::value(&QLMSpin), "value of the spin sites for the QLM [default 1/2]")
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
         ("H_Jxm"  , "nearest-neighbor hopping in -x-direction")
         ("H_Jxp"  , "nearest-neighbor hopping in +x-direction")
         ("H_Jym"  , "nearest-neighbor hopping in -y-direction")
         ("H_Jyp"  , "nearest-neighbor hopping in +y-direction")
         ("H_Jx"   , "nearest-neighbor hopping in x-direction")
         ("H_Jy"   , "nearest-neighbor hopping in y-direction")
         ("H_J"    , "nearest-neighbor hopping")
         ("H_U"    , "on-site Coulomb interaction N*(N-1)/2")
         ("H_Vx"   , "nearest-neighbour Coulomb repulsion in x-direction")
         ("H_Vy"   , "nearest-neighbour Coulomb repulsion in x-direction")
         ("H_V"    , "nearest-neighbour Coulomb repulsion")
         ("H_delta", "QLM gauge site potential")
         ("H_alpha", "QLM matter site interaction")
         ("H_W"    , "QLM gauge site repulsion")
         ("H_chi"  , "QLM staggering potential")
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
      UnitCell Cell(repeat(Site, 3*CellSize));
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2"), I(Cell, "I");

      // Define hopping terms and near-neighbour interactions.
      UnitCellMPO Jxm, Jxp, Jym, Jyp, Vx, Vy;

      // the XY configuration is special
      if (x == 0)
      {
         for (int i = 0; i < y; ++i)
         {
            Jxm += BH(0)[3*i]*B(0)[3*i+1];
            Jxp += B(0)[3*i]*BH(0)[3*i+1];
            Jxm += BH(0)[3*i+1]*B(1)[3*i];
            Jxp += B(0)[3*i+1]*BH(1)[3*i];
            Jym += BH(0)[3*i]*B(0)[3*i+2];
            Jyp += B(0)[3*i]*BH(0)[3*i+2];
            Jym += BH(0)[3*i+2]*B(0)[(3*i+3)%(3*y)];
            Jyp += B(0)[3*i+2]*BH(0)[(3*i+3)%(3*y)];
            Vx += N(0)[3*i]*N(0)[3*i+1];
            Vx += N(0)[3*i+1]*N(1)[3*i];
            Vy += N(0)[3*i]*N(0)[3*i+2];
            Vy += N(0)[3*i+2]*N(0)[(3*i+3)%(3*y)];
         }
      }
      else
      {
         // ...
      }

      Lattice["H_Jxm"] = sum_unit(Jxm);
      Lattice["H_Jxp"] = sum_unit(Jxp);
      Lattice["H_Jx"] = sum_unit(Jxm+Jxp);
      Lattice.func("H_Jx")(arg("theta")) = "exp(-i*theta)*H_Jxm + exp(i*theta)*H_Jxp";
      Lattice["H_Jym"] = sum_unit(Jym);
      Lattice["H_Jyp"] = sum_unit(Jyp);
      Lattice["H_Jy"] = sum_unit(Jym+Jyp);
      Lattice.func("H_Jy")(arg("theta")) = "exp(-i*theta)*H_Jym + exp(i*theta)*H_Jyp";
      Lattice["H_J"] = sum_unit(Jxp+Jxm+Jyp+Jym);

      Lattice["H_Vx"] = sum_unit(Vx);
      Lattice["H_Vy"] = sum_unit(Vy);
      Lattice["H_V"] = sum_unit(Vx+Vy);

      // Define on-site interactions.
      UnitCellMPO U;

      for (int i = 0; i < 3*y; ++i)
         U += 0.5*N2(0)[i];

      Lattice["H_U"] = sum_unit(U);

      if (x == 0)
      {
         // Define QLM potentials.
         UnitCellMPO delta, alpha, W;

         for (int i = 0; i < y; ++i)
         {
            delta += N(0)[3*i+1] + N(0)[3*i+2];
            alpha += 0.5*N2(0)[3*i];
            W += N(0)[3*i+1] * N(0)[3*i+2];
            W += N(0)[3*i+1] * N(0)[(3*i-1+3*y)%(3*y)];
            W += N(0)[3*i+1] * N(-1)[3*i+1];
            W += N(0)[3*i+2] * N(0)[(3*i-1+3*y)%(3*y)];
            W += N(0)[3*i+2] * N(-1)[3*i+1];
            W += N(0)[(3*i-1+3*y)%(3*y)] * N(-1)[3*i+1];
         }

         Lattice["H_delta"] = sum_unit(delta, 3*y);
         Lattice["H_alpha"] = sum_unit(alpha, 3*y);
         Lattice["H_W"] = sum_unit(W, 3*y);

         if (y%2 == 0)
         {
            UnitCellMPO chi;
            for (int i = 0; i < y; ++i)
            {
               chi += N(0)[3*i+1] + N(0)[3*i+2];
               chi -= N(1)[3*i+1] + N(1)[3*i+2];
               ++i;
               chi -= N(0)[3*i+1] + N(0)[3*i+2];
               chi -= N(1)[3*i+1] + N(1)[3*i+2];
            }
            Lattice["H_chi"] = sum_unit(0.5*chi, 6*y);
         }

         // Define Gauss's law operators.
         // Note: Ideally, this operator should be multiplied by -1 for an antimatter site.
         UnitCellMPO G[y];

         for (int i = 0; i < y; ++i)
         {
            G[i] = N(0)[3*i] + 0.5 * (N(0)[3*i+1] + N(0)[3*i+2] + N(0)[(3*i-1+3*y)%(3*y)] + N(-1)[3*i+1]) - 4.0 * QLMSpin * I(0)[3*i];
            G[i].set_description("Gauss's law operator for matter site " + std::to_string(3*i));
            Lattice.GetUnitCell().assign_operator("G" + std::to_string(3*i), G[i]);
         }
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
