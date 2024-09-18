// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// models/hubbardcylinder-u1su2-k.cpp
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
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
#include "models/fermion-u1su2-k.h"
#include "common/terminal.h"
#include "common/prog_options.h"

namespace prog_opt = boost::program_options;

int const NN = 8;

int main(int argc, char** argv)
{
   try
   {
      std::string FileName;
      int w = NN;

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
      OpDescriptions.set_description("U(1)xSU(2) Fermi Hubbard model");
      OpDescriptions.author("IP McCulloch", "ianmcc@physics.uq.edu.au");
      OpDescriptions.add_operators()
         ("H_tx" , "nearest neighbor hopping in y-direction")
         ("H_ty" , "nearest neighbor hopping in x-direction")
         ("H_t"  , "nearest neighbor hopping")
         ("H_U"  , "on-site Coulomb interaction n_up*n_down")
         ("H_Us" , "on-site Coulomb interaction (n_up-1/2)(n_down-1/2)")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Operators:\n" << OpDescriptions;
         return 1;
      }

      std::vector<LatticeSite> Sites;
      for (int k = 0; k < w; ++k)
      {
         Sites.push_back(FermionU1SU2_K<NN>(k));
      }

      std::vector<LatticeSite> Sites2 = {Sites[0], Sites[4], Sites[1], Sites[5], Sites[2], Sites[6],
                                         Sites[3], Sites[7]};

      //std::vector<int> kk = {0,2,4,6,1,3,5,7};
      std::vector<int> kk = {0,1,2,3,4,5,6,7};

      UnitCell Cell(Sites);
      InfiniteLattice Lattice(&Cell);
      UnitCellOperator CH(Cell, "CH"), C(Cell, "C"), Pdouble(Cell, "Pdouble"),
         Hu(Cell, "Hu"), N(Cell, "N");

      UnitCellMPO tx, ty, Pd, hu;

      QuantumNumbers::QNConstructor<QuantumNumbers::U1,QuantumNumbers::SU2,
         QuantumNumbers::Zn<NN>> QN(Cell.GetSymmetryList());

      for (int k = 0; k < w; ++k)
      {
         ty += std::cos(math_const::pi*2*k / w) * N(0)[k];
      }

      for (int k = 0; k < w; ++k)
      {
         tx += dot(CH(0)[k], C(1)[k]) + dot(C(0)[k], CH(1)[k]);
      }

      for (int k = 0; k < w; ++k)
      {
         for (int l = 0; l < w; ++l)
         {
            for (int m = 0; m < w; ++m)
            {
               UnitCellMPO Op = (1.0/w) * prod(CH(0)[kk[k]], C(0)[kk[l]], QN(0,0,k-l))
                  * prod(CH(0)[kk[m]], C(0)[kk[(k-l+m+w)%w]], QN(0,0,l-k));
               CHECK(Op.size() != 0);
               hu += Op;
            }
         }
         Pd -= N(0)[k];
      }
      Pd += hu;

      Lattice["H_tx"] = sum_unit(tx);
      Lattice["H_ty"] = sum_unit(ty);
      Lattice["H_t"] = sum_unit(tx+ty);
      Lattice["H_Us"] = sum_unit(Pd);
      Lattice["H_U"] = sum_unit(hu);

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
