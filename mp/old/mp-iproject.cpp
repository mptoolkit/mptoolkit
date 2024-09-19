// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-iproject.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

#include "mps/infinitewavefunction.h"
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

// construct a projector onto a single quantum number subspace
MatrixOperator
ProjectOnto(VectorBasis const& b, QuantumNumbers::QuantumNumber const& q)
{
   // firstly, find out the dimension of the q subspace in b
   int Dim = 0;
   unsigned s = 0;
   while (s < b.size() && b[s] != q)
   {
      TRACE(b[s]);
      ++s;
   }
   if (s == b.size())
   {
      std::cerr << "fatal: quantum number doesn't exist in the basis, wavefunction is zero!\n";
      exit(1);
   }

   Dim = b.dim(s);

   VectorBasis New(b.GetSymmetryList());
   New.push_back(q, Dim);

   MatrixOperator Result(New, b);
   Result(0, s) = LinearAlgebra::identity_matrix<double>(Dim);
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      std::string q;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("input-wavefunction", prog_opt::value<std::string>(), "input wavefunction (required)")
         ("quantum-number", prog_opt::value(&q), "input wavefunction (required)")
         ;

      prog_opt::positional_options_description p;
      p.add("input-wavefunction", 1);
      p.add("quantum-number", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("quantum-number") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-project input-wavefunction quantum-number\n";
         std::cerr << desc << "\n";
         return 1;
      }

      std::string Wavefunc = vm["input-wavefunction"].as<std::string>();

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(Wavefunc, mp_pheap::CacheSize());

      InfiniteWavefunction Psi = *PsiPtr.mutate();

      // construct projector onto the quantum number

      QuantumNumbers::QuantumNumber Q(Psi.Psi.GetSymmetryList(), q);

      TRACE(Q);

      MatrixOperator Projector = ProjectOnto(Psi.Psi.Basis2(), Q);

      Psi.C_right = Projector * Psi.C_right * herm(Projector);
      Psi.Psi.set_back(prod(Psi.Psi.get_back(), herm(Projector)));

      MatrixOperator ShiftedProjector = delta_shift(Projector, Psi.QShift);

      Psi.C_old = ShiftedProjector * Psi.C_old * herm(ShiftedProjector);
      Psi.Psi.set_front(prod(ShiftedProjector, Psi.Psi.get_front()));

      orthogonalize(Psi);

      PsiPtr = new InfiniteWavefunction(Psi);

      pheap::ShutdownPersistent(PsiPtr);
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
