// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-scale.cpp
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

#include "wavefunction/mpwavefunction.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "parser/number-parser.h"

#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"

namespace prog_opt = boost::program_options;
using formatting::format_complex;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string PsiStr;
      std::string Expression;
      int UnitCell = 0;
      bool Log = false;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("unitcell,u", prog_opt::value(&UnitCell), "For iMPS, interpret the scaling factor as per this unit cell size")
         ("log", prog_opt::bool_switch(&Log), "For an iMPS, interpret the value as the logarithm of the scale factor")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity (can be used more than once)")
         ("help", "show this help message")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value<std::string>(&PsiStr), "psi")
         ("expr", prog_opt::value<std::string>(&Expression), "expr")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("expr", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("expr") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi1> <factor>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize());

      std::complex<double> x = ParseNumber(Expression);

      if (Psi->is<FiniteWavefunctionLeft>())
      {
         if (vm.count("log"))
         {
            x = std::exp(x);
         }
         if (vm.count("unitcell"))
         {
            std::cerr << basename(argv[0]) << ": error: <unitcell> option is not valid for finite wavefunctions.\n";
            return 1;
         }
         Psi.mutate()->get<FiniteWavefunctionLeft>() *= x;
      }
      else if (Psi->is<InfiniteWavefunctionLeft>())
      {
         if (vm.count("log"))
         {
            if (vm.count("unitcell"))
            {
               int u = Psi->get<InfiniteWavefunctionLeft>().size();
               x = x * (double(u) / double(UnitCell));
            }
            Psi.mutate()->get<InfiniteWavefunctionLeft>().scale_log(x);
         }
         else
         {
            if (vm.count("unitcell"))
            {
               int u = Psi->get<InfiniteWavefunctionLeft>().size();
               x = std::pow(x, double(u) / double(UnitCell));
               if (Verbose > 0)
               {
                  std::cout << "Rescaling to unit cell " << UnitCell << " with scaling factor "
                     << format_complex(x) << '\n';
               }
            }
            Psi.mutate()->get<InfiniteWavefunctionLeft>() *= x;
         }
      }
      else if (Psi->is<InfiniteWavefunctionRight>())
      {
         if (vm.count("log"))
         {
            if (vm.count("unitcell"))
            {
               int u = Psi->get<InfiniteWavefunctionRight>().size();
               x = x * (double(u) / double(UnitCell));
            }
            Psi.mutate()->get<InfiniteWavefunctionRight>().scale_log(x);
         }
         else
         {
            if (vm.count("unitcell"))
            {
               int u = Psi->get<InfiniteWavefunctionRight>().size();
               x = std::pow(x, double(u) / double(UnitCell));
               if (Verbose > 0)
               {
                  std::cout << "Rescaling to unit cell " << UnitCell << " with scaling factor "
                     << format_complex(x) << '\n';
               }
            }
            Psi.mutate()->get<InfiniteWavefunctionRight>() *= x;
         }
      }
      else
      {
         std::cerr << "mp-scale: fatal: wavefunction type " << Psi->Type() << " is not supported.\n";
         return 1;
      }
      Psi.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));

      pheap::ShutdownPersistent(Psi);
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
