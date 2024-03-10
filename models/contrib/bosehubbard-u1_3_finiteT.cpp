// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/contrib/bosehubbard-u1_3_finiteT.cpp
//
// Copyright (C) 2015-2021 Ian McCulloch <ian@qusim.net>
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
#include "models/contrib/boson-u1-aux.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      int MaxN = DefaultMaxN;
      int L = 0;
      double Beta = 1.0;
      double Sigma = 0.02;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&MaxN),
          FormatDefault("Maximum number of bosons per site", MaxN).c_str())
         ("Unit cell length,l", prog_opt::value(&L),
          FormatDefault("Unit cell length", L).c_str())
         ("beta", prog_opt::value(&Beta),
          FormatDefault("Trap height", Beta).c_str())
         ("sigma", prog_opt::value(&Sigma),
          FormatDefault("Trap width", Sigma).c_str())
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.set_description("U(1) Bose-Hubbard model");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_J"			, "nearest-neighbor hopping\n")
         ("H_U"			, "on-site Coulomb repulsion N*(N-1)/2\n")
         ("H_V"			, "nearest-neighbour Coulomb repulsion\n")
         ("H_N"			, "on-site potential\n")
         ("H_trap"		, "Trap potential\n")
         ("H_trap_int"		, "Trap potential - interaction part\n")
         ("H_trap_int2"		, "Trap potential - interaction part - TG limit\n")
	 ("H_harmonic"		, "Harmonic trap potential\n")
         ("H_J2"		, "nearest-neighbor hopping for auxiliary sites\n")
	 ("IP"			, "Artificial projector Hamiltonian\n")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << OpDescriptions << '\n';
         return 1;
      }

      LatticeSite RealSite = BosonU1(MaxN, "NReal"); // Real sites
      LatticeSite AuxSite = BosonU1Aux(MaxN, "NAux"); // Auxiliary sites
      UnitCell Cell(SymmetryList("NReal:U(1),NAux:U(1)"), RealSite,AuxSite);
      UnitCellOperator BH(Cell, "BH"), B(Cell, "B"), N(Cell, "N"), N2(Cell, "N2");
      InfiniteLattice Lattice(&Cell);

      UnitCellMPO H_t, H_t_int, H_t_int2, H_t_x2;
      for (int i = 0; i < L; i++)
//      for (int i = 0; i <= (L/2)-1; i++)
      {
	double Ramp = (i - (L+1.0)/2.0) / (L-1.0);
//	double Ramp = (i - ((L/2.0)+1.0)/2.0) / ((L/2.0)-1.0);
	H_t += ((Beta/(2*pow(Sigma,4))) * (Ramp*Ramp - Sigma*Sigma) / (Beta + exp(Ramp*Ramp/(2*Sigma*Sigma)))) * N(i)[0];
	H_t_int += -pow((1 + Beta*exp(-Ramp*Ramp / (2*Sigma*Sigma))), 2) * N(i)[0];
	H_t_int2 += -pow((1 + Beta*exp(-Ramp*Ramp / (2*Sigma*Sigma))), 4) * N(i)[0];
	H_t_x2 += Ramp*Ramp * N(i)[0];
      }

      // Hamiltonian for real sites
      Lattice["H_J"] = sum_unit(BH(0)[0]*B(1)[0] + B(0)[0]*BH(1)[0]);
      Lattice["H_U"] = sum_unit(0.5*N2(0)[0]);
      Lattice["H_V"] = sum_unit(N(0)[0]*N(1)[0]);
      Lattice["H_N"] = sum_unit(N(0)[0]);
//      Lattice["H_trap"] = sum_unit(H_t, L);
//      Lattice["H_trap_int"] = sum_unit(H_t_int, L);
//      Lattice["H_trap_int2"] = sum_unit(H_t_int2, L);
//      Lattice["H_harmonic"] = sum_unit(H_t_x2, L);
      Lattice["H_trap"] = sum_unit(H_t, 2*L); // Use 2*L instead of L becuase of the coarse-graining of the real and auxiliary sites during real- and imaginary-time evolution
      Lattice["H_trap_int"] = sum_unit(H_t_int, 2*L); // Use 2*L instead of L becuase of the coarse-graining of the real and auxiliary sites during real- and imaginary-time evolution
      Lattice["H_trap_int2"] = sum_unit(H_t_int2, 2*L); // Use 2*L instead of L becuase of the coarse-graining of the real and auxiliary sites during real- and imaginary-time evolution
      Lattice["H_harmonic"] = sum_unit(H_t_x2, 2*L); // Use 2*L instead of L becuase of the coarse-graining of the real and auxiliary sites during real- and imaginary-time evolution

      // Hamiltonian for auxiliary sites
      Lattice["H_J2"] = sum_unit(BH(0)[1]*B(1)[1] + B(0)[1]*BH(1)[1]);

      UnitCellOperator P_0(Cell, "P_0"), P_1(Cell, "P_1"), P_2(Cell, "P_2"), P_3(Cell, "P_3"), P_4(Cell, "P_4"),
		       P_5(Cell, "P_5"), P_6(Cell, "P_6"), P_7(Cell, "P_7"), P_8(Cell, "P_8"), P_9(Cell, "P_9"),
		       P_10(Cell, "P_10");
/*
      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1] + P_4(0)[0]*P_4(0)[1] + P_5(0)[0]*P_5(0)[1]
			+ P_6(0)[0]*P_6(0)[1] + P_7(0)[0]*P_7(0)[1] + P_8(0)[0]*P_8(0)[1]
			+ P_9(0)[0]*P_9(0)[1]);
*/
/*
      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1] + P_4(0)[0]*P_4(0)[1] + P_5(0)[0]*P_5(0)[1]
			+ P_6(0)[0]*P_6(0)[1] + P_7(0)[0]*P_7(0)[1] + P_8(0)[0]*P_8(0)[1]);
*/
/*
      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1] + P_4(0)[0]*P_4(0)[1] + P_5(0)[0]*P_5(0)[1]
			+ P_6(0)[0]*P_6(0)[1] + P_7(0)[0]*P_7(0)[1]);
*/
/*
      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1] + P_4(0)[0]*P_4(0)[1] + P_5(0)[0]*P_5(0)[1]
			+ P_6(0)[0]*P_6(0)[1]);
*/

      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1] + P_4(0)[0]*P_4(0)[1] + P_5(0)[0]*P_5(0)[1]);

/*
      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1] + P_4(0)[0]*P_4(0)[1]);
*/
/*
      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]
			+ P_3(0)[0]*P_3(0)[1]);
*/
//      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1] + P_2(0)[0]*P_2(0)[1]);

//      Lattice["IP"] = (1/(MaxN+1))*sum_unit(P_0(0)[0]*P_0(0)[1] + P_1(0)[0]*P_1(0)[1]);


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
