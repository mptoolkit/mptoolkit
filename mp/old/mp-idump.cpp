// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-idump.cpp
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

#include "matrixproduct/infinitewavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

using namespace LinearAlgebra;

void DumpMatrixOperator(MatrixOperator const& M, bool ShowReal, bool ShowImag)
{
   for (const_iterator<MatrixOperator>::type I = iterate(M); I; ++I)
   {
      for (const_inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
      {
         if (M.Basis1().size() != 1 || M.Basis2().size() != 1)
            std::cout << "Component at (" << J.index1() << ", " << J.index2() << ")\n";
         if (!ShowReal && !ShowImag)
            std::cout << *J << '\n';
         if (ShowReal)
            std::cout << real(*J) << '\n';
         if (ShowImag)
            std::cout << imag(*J) << '\n';
      }
   }
}

void DumpMPState(MPStateComponent const& A, bool ShowReal, bool ShowImag)
{
   for (unsigned i = 0; i < A.size(); ++i)
   {
      std::cout << "Component A[" << i << "], with quantum number " << A.LocalBasis()[i] << '\n';
      DumpMatrixOperator(A[i], ShowReal, ShowImag);
   }
}

void DumpWavefunction(LinearWavefunction const& Psi, bool ShowReal, bool ShowImag)
{
   LinearWavefunction::const_iterator I = Psi.begin();
   int i = 1;
   while (I != Psi.end())
   {
      if (Psi.size() != 1)
         std::cout << "Wavefunction at position " << i << " in the unit cell:\n";
      DumpMPState(*I, ShowReal, ShowImag);
      ++I; ++i;
   }
}

void DumpDiagonalDense(Matrix<std::complex<double> > const& M, bool ShowReal, bool ShowImag)
{
   CHECK_EQUAL(size1(M), size2(M));
   for (unsigned i = 0; i < size1(M); ++i)
   {
      if (!ShowReal && !ShowImag)
         std::cout << M(i,i);
      else
      {
         if (ShowReal)
            std::cout << real(M(i,i));
         if (ShowImag)
            std::cout << imag(M(i,i));
      }
      std::cout << '\n';
   }
}

void DumpDiagonalMatrix(MatrixOperator const& M, bool ShowReal, bool ShowImag)
{
   for (unsigned i = 0; i < M.Basis1().size(); ++i)
   {
      const_inner_iterator<MatrixOperator>::type I = iterate_at(M.data(), i, i);
      if (I)
      {
         if (M.Basis1().size() != 1)
            std::cout << "Component at (" << i << ", " << i << ")\n";
         DumpDiagonalDense(*I, ShowReal, ShowImag);
      }
   }
}

int main(int argc, char** argv)
{
   try
   {
      bool ShowReal = false;
      bool ShowImag = false;
      std::string PsiStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("real,r", prog_opt::bool_switch(&ShowReal),
          "display only the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImag),
          "display only the imaginary part of the result")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&PsiStr), "psi")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("lhs") == 0 || (ShowReal && ShowImag))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-idump [options] <psi>\n";
         std::cerr << desc << '\n';
         if (ShowReal && ShowImag)
            std::cerr << "--real and --imag options are mutually exclusive.\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      pvalue_ptr<InfiniteWavefunction> Psi
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      DumpWavefunction(Psi->Psi, ShowReal, ShowImag);
      std::cout << "\nSingular value matrix Lambda:\n";
      DumpDiagonalMatrix(Psi->C_right, ShowReal, ShowImag);

      pheap::Shutdown();
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
