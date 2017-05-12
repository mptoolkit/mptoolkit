// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-overlap.cpp
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

namespace prog_opt = boost::program_options;

void PrintFormat(std::complex<double> Value, bool ShowRealPart, bool ShowImagPart,
                 bool ShowMagnitude, bool ShowArgument, bool ShowRadians)
{
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(Value) << "    ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
         Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
   std::cout << std::endl;
}


int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      bool UseTempFile = false;
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false, ShowCorrLength = false;
      std::string LhsStr, RhsStr;
      bool Quiet = false;
      bool Reflect = false;
      bool Conj = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("cart,c", prog_opt::bool_switch(&ShowCartesian),
          "show the result in cartesian coordinates [equivalent to --real --imag]")
         ("polar,p", prog_opt::bool_switch(&ShowPolar),
          "show the result in polar coodinates [equivalent to --mag --arg]")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "display the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "display the argument in radians instead of degrees")
         ("reflect", prog_opt::bool_switch(&Reflect),
          "reflect psi2")
         ("conj", prog_opt::bool_switch(&Conj),
          "complex conjugate psi2")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the column headings")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
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
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi1> <psi2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (!Quiet)
         print_preamble(std::cout, argc, argv);

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
          && !ShowCartesian && !ShowPolar && !ShowArgument
          && !ShowRadians && !ShowCorrLength)
      {
         ShowCartesian = true;
         ShowPolar = true;
         ShowCorrLength = true;
      }

      if (ShowCartesian)
      {
         ShowRealPart = true;
         ShowImagPart = true;
      }
      if (ShowPolar)
      {
         ShowMagnitude = true;
         ShowArgument = true;
      }
      if (ShowRadians)
         ShowArgument = true;

      if (Verbose)
         std::cout << "Loading LHS wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi1Ptr;
      if (UseTempFile)
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1Ptr = pheap::ImportHeap(LhsStr);
      }
      else
      {
         Psi1Ptr = pheap::OpenPersistent(LhsStr, mp_pheap::CacheSize(), true);
      }
      FiniteWavefunctionLeft Psi1 = Psi1Ptr->get<FiniteWavefunctionLeft>();

      if (Verbose)
         std::cout << "Loading RHS wavefunction...\n";

      int Size = Psi1.size();

      pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(RhsStr);
      FiniteWavefunctionLeft Psi2 = Psi2Ptr->get<FiniteWavefunctionLeft>();

      if (ExtractLocalBasis(Psi1) != ExtractLocalBasis(Psi2))
      {
         std::cerr << "mp-ioverlap: fatal: local basis for LHS and RHS wavefunctions do not match!\n";
         return 1;
      }

      if (Verbose)
         std::cout << "Calculating overlap...\n";

      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi2..." << std::endl;
         inplace_reflect(Psi2);
      }

      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi2..." << std::endl;
         inplace_conj(Psi2);
      }


      if (!Quiet)
      {
         std::cout << "#size of wavefunction is " << Size << " sites\n";
         if (ShowRealPart)
            std::cout << "#real                   ";
         if (ShowImagPart)
            std::cout << "#imag                   ";
         if (ShowMagnitude)
            std::cout << "#magnitude              ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
         std::cout << '\n';
         std::cout << std::left;
      }

      std::complex<double> x = overlap(Psi1, Psi2);

      PrintFormat(x, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);

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
