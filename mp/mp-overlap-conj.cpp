// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-overlap-conj.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      bool ShowRealPart = false, ShowImagPart = false;
      std::string LhsStr, RhsStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("r,real", prog_opt::bool_switch(&ShowRealPart),
          "display only the real part of the result")
         ("i,imag", prog_opt::bool_switch(&ShowImagPart),
          "display only the imaginary part of the result")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&LhsStr), "psi1")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("rhs") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-overlap-conj [options] <psi1> <psi2>\n";
         std::cerr << "Calculates the inner product <psi1|conj(psi2)>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<MPWavefunction> Psi1 = pheap::OpenPersistent(LhsStr, CacheSize, true);
      pvalue_ptr<MPWavefunction> Psi2 = pheap::ImportHeap(RhsStr);

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::complex<double> x = overlap_conj(*Psi1, *Psi2);

      if (ShowRealPart || ShowImagPart)
      {
         if (ShowRealPart)
         {
            std::cout << x.real();
            if (ShowImagPart)
               std::cout << "   " << x.imag();
         }
         else // if we get here then ShowImagPart is true and ShowRealPart is false
            std::cout << x.imag();
      }
      else // default to C++ complex output
         std::cout << x;

      std::cout << '\n';

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
