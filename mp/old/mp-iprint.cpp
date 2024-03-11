// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-iprint.cpp
//
// Copyright (C) 2014-2017 Ian McCulloch <ian@qusim.net>
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

#include "wavefunction/mpwavefunction.h"
#include "tensor/tensor_eigen.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("input-wavefunction", prog_opt::value<std::string>(), "input wavefunction (required)")
         ;

      prog_opt::positional_options_description p;
      p.add("input-wavefunction", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("input-wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-iprint [options] input-wavefunction\n";
         std::cerr << desc << "\n";
         return 1;
      }

      std::string Wavefunc = vm["input-wavefunction"].as<std::string>();

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(Wavefunc, mp_pheap::CacheSize(), true);

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

#if 0

      std::cout << "Symmetry list=" << Psi->C_right.GetSymmetryList() << '\n';
      std::cout << "Transforms as=" << Psi->shift() << '\n';

      std::cout << "Number of states=" << Psi->C_right.Basis1().total_dimension() << '\n';
      std::cout << "Degree=" << Psi->C_right.Basis1().total_degree() << '\n';
      std::cout << "Unit cell size=" << Psi->Psi.size() << '\n';

      std::cout << "Orthogonality fidelity=" << (1.0 - orthogonality_fidelity(*Psi)) << '\n';

      InfiniteWavefunction P = *Psi;

      std::cout << "\nC_old = " << P.C_old << '\n';
      std::cout << "\nC_right = " << P.C_right << '\n';

      bool Symmetric = true;
      if (Symmetric)
      {
         MatrixOperator LambdaSqrt = SqrtDiagonal(P.C_old);
         MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);
         P.Psi.set_back(prod(P.Psi.get_back(), delta_shift(LambdaSqrt, adjoint(P.QShift))));
         P.Psi.set_front(prod(LambdaInvSqrt, P.Psi.get_front()));
      }
#endif

      int Site = 0;
      for (auto const& x : Psi)
      {
         std::cout << "\n***** Site " << (Site++) << " *****\n\n";
         std::cout << x << "\n";
      }

   }
   catch(std::exception& e)
   {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
   }
   catch(...)
   {
      std::cerr << "Exception of unknown type!\n";
   }
   return 0;
}
