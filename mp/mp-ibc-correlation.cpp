// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-correlation.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "lattice/infinite-parser.h"

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
      bool Quiet = false;
      bool Conj = false;
      //bool Simple = false;
      std::string String;
      ProductMPO StringOp;

      std::string InputPrefix;
      std::string TimestepStr;
      std::string BetastepStr;
      std::string InitialTimeStr;
      std::string InitialBetaStr;
      int N = 0;
      int XMin = 0;
      int XMax = 0;
      int YMax = 0;
      int UCSize = -1;
      int Digits = 0;
      std::complex<double> InitialTime = 0.0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("cart,c", prog_opt::bool_switch(&ShowCartesian),
          "Show the result in cartesian coordinates [equivalent to --real --imag]")
         ("polar,p", prog_opt::bool_switch(&ShowPolar),
          "Show the result in polar coodinates [equivalent to --mag --arg]")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "Display the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "Display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "Display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "Display the argument (angle) of the result")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "Display the argument in radians instead of degrees")
         ("conj", prog_opt::bool_switch(&Conj),
          "Complex conjugate psi2")
         //("simple", prog_opt::bool_switch(&Simple),
         // "Assume that psi1 and psi2 have the same semi-infinite boundaries")
         ("string", prog_opt::value(&String),
          "Use this string MPO representation for the cylinder translation operator")
	 ("wavefunction,w", prog_opt::value(&InputPrefix), "Prefix for input wavefunctions")
	 ("timestep,t", prog_opt::value(&TimestepStr), "Timestep")
	 ("precision", prog_opt::value(&Digits), "Decimal precision in time value of the filenames")
	 ("num-timesteps,n", prog_opt::value(&N), FormatDefault("Number of timesteps", N).c_str())
	 ("xmin", prog_opt::value(&XMin), FormatDefault("Minimum value of x", XMin).c_str())
	 ("xmax,x", prog_opt::value(&XMax), FormatDefault("Maximum value of x", XMax).c_str())
	 ("ymax,y", prog_opt::value(&YMax), FormatDefault("Maximum value of y (requires --string)", YMax).c_str())
	 ("ucsize", prog_opt::value(&UCSize), "Unit cell size [default left boundary unit cell size]")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "Don't show the column headings")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "Extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);


      if (vm.count("help") > 0 || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]" << std::endl;
         std::cerr << desc << std::endl;
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

      std::complex<double> Timestep = 0.0;

      if (!TimestepStr.empty())
         Timestep += ParseNumber(TimestepStr);

      Digits = std::max(Digits, formatting::digits(Timestep));

      if (Verbose)
         std::cout << "Loading RHS wavefunction..." << std::endl;

      std::string InitialFilename = InputPrefix;

      if (std::real(Timestep) != 0.0)
         InitialFilename += ".t" + formatting::format_digits(InitialTime, Digits);

      if (std::imag(Timestep) != 0.0)
         InitialFilename += ".b" + formatting::format_digits(InitialTime, Digits);

      pvalue_ptr<MPWavefunction> Psi1Ptr;
      if (UseTempFile)
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1Ptr = pheap::ImportHeap(InitialFilename);
      }
      else
      {
         Psi1Ptr = pheap::OpenPersistent(InitialFilename, mp_pheap::CacheSize(), true);
      }
      IBCWavefunction Psi1 = Psi1Ptr->get<IBCWavefunction>();

      // If UCSize is not specified, then we set it to the left boundary unit
      // cell size.
      if (UCSize = -1)
         UCSize = Psi1.Left.size();
      
      if (vm.count("string"))
      {
         InfiniteLattice Lattice;
         std::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(String);
      }
      else
      {
         StringOp = ProductMPO::make_identity(ExtractLocalBasis(Psi1.Left));
      }

      // Get the string operators representing the rotations by 0, 1, ..., YMax sites.
      std::vector<ProductMPO> StringOpVec;
      StringOpVec.push_back(ProductMPO::make_identity(ExtractLocalBasis(Psi1.Left)));

      for (int y = 1; y <= YMax; ++y)
         StringOpVec.push_back(pow(StringOp, y));

      IBCWavefunction Psi2;

      Psi2 = Psi1;
      if (Conj)
         inplace_conj(Psi1);

      std::vector<StateComponent> EVec, FVec;
      // Scaling factor to remove the spurious phase when using the complex conjugate.
      std::vector<std::complex<double>> PhaseFactorVec;
      
      auto Op = StringOpVec.begin();
      for (int y = 0; y <= YMax; ++y)
      {
         StateComponent E, F;
         std::tie(E, F) = get_boundary_transfer_eigenvectors(Psi2, *Op, Psi1, Verbose);
         EVec.push_back(E);
         FVec.push_back(F);

         std::complex<double> Overlap = overlap(Psi2, Psi1, E, F);
         // TODO: Here we assume that the phases for the x = t = 0 overlaps
         // should all be 0: this may not be the case, and we need to get the
         // correct phases from somewhere.
         PhaseFactorVec.push_back(std::abs(Overlap) / Overlap);

         ++Op;
      }

      if (!Quiet)
      {
         if (std::real(Timestep) != 0.0)
            std::cout << "#time         ";
         if (std::imag(Timestep) != 0.0)
            std::cout << "#beta         ";
         if (XMin != 0 || XMax != 0)
            std::cout << "#x            ";
         if (YMax != 0)
            std::cout << "#y            ";
         if (ShowRealPart)
            std::cout << "#real                   ";
         if (ShowImagPart)
            std::cout << "#imag                   ";
         if (ShowMagnitude)
            std::cout << "#magnitude              ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
         std::cout << std::endl;
         std::cout << std::left;
      }

      for (int tstep = 0; tstep <= N; ++tstep)
      {
         if (Verbose)
            std::cout << "Loading LHS wavefunction..." << std::endl;

         std::string Filename = InputPrefix;
         std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), Digits);
         std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), Digits);

         if (std::real(Timestep) != 0.0)
            Filename += ".t" + TimeStr;

         if (std::imag(Timestep) != 0.0)
            Filename += ".b" + BetaStr;

         pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(Filename);
         Psi2 = Psi2Ptr->get<IBCWavefunction>();

         Psi2.WindowOffset += (XMin-1) * UCSize;

         for (int x = XMin; x <= XMax; ++x)
         {
            Psi2.WindowOffset += UCSize;

            Op = StringOpVec.begin();
            auto E = EVec.begin();
            auto F = FVec.begin();
            auto PhaseFactor = PhaseFactorVec.begin();

            for (int y = 0; y <= YMax; ++y)
            {
               std::complex<double> Overlap;

               Overlap = *PhaseFactor * overlap(Psi2, *Op, Psi1, *E, *F, Verbose);

               if (std::real(Timestep) != 0.0)
                  std::cout << std::setw(10) << TimeStr << "    ";
               if (std::imag(Timestep) != 0.0)
                  std::cout << std::setw(10) << BetaStr << "    ";
               if (XMin != 0 || XMax != 0)
                  std::cout << std::setw(10) << x << "    ";
               if (YMax != 0)
                  std::cout << std::setw(10) << y << "    ";

               PrintFormat(Overlap, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);

               ++Op, ++E, ++F, ++PhaseFactor;
            }
         }
      }
      
      if (Conj)
      {
         std::complex<double> FinalTime = double(N) * Timestep;
            Psi2.WindowOffset -= XMax * UCSize;

         for (int tstep = 1; tstep <= N; ++tstep)
         {
            if (Verbose)
               std::cout << "Loading RHS wavefunction..." << std::endl;

            // The time strings for the input filename.
            std::string Filename = InputPrefix;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), Digits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), Digits);

            if (std::real(Timestep) != 0.0)
               Filename += ".t" + TimeStr;

            if (std::imag(Timestep) != 0.0)
               Filename += ".b" + BetaStr;

            pvalue_ptr<MPWavefunction> Psi1Ptr = pheap::ImportHeap(Filename);
            IBCWavefunction Psi1 = Psi1Ptr->get<IBCWavefunction>();
            inplace_conj(Psi1);

            // The time strings for outputting to the command line.
            TimeStr = formatting::format_digits(std::real(InitialTime + FinalTime + double(tstep)*Timestep), Digits);
            BetaStr = formatting::format_digits(-std::imag(InitialTime + FinalTime + double(tstep)*Timestep), Digits);

            Psi2.WindowOffset += (XMin-1) * UCSize;

            for (int x = XMin; x <= XMax; ++x)
            {
               Psi2.WindowOffset += UCSize;

               Op = StringOpVec.begin();
               auto E = EVec.begin();
               auto F = FVec.begin();
               auto PhaseFactor = PhaseFactorVec.begin();

               for (int y = 0; y <= YMax; ++y)
               {
                  std::complex<double> Overlap;

                  Overlap = *PhaseFactor * overlap(Psi2, *Op, Psi1, *E, *F, Verbose);

                  if (std::real(Timestep) != 0.0)
                     std::cout << std::setw(10) << TimeStr << "    ";
                  if (std::imag(Timestep) != 0.0)
                     std::cout << std::setw(10) << BetaStr << "    ";
                  if (XMin != 0 || XMax != 0)
                     std::cout << std::setw(10) << x << "    ";
                  if (YMax != 0)
                     std::cout << std::setw(10) << y << "    ";

                  PrintFormat(Overlap, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);

                  ++Op, ++E, ++F, ++PhaseFactor;
               }
            }
            Psi2.WindowOffset -= XMax * UCSize;
         }
      }

      pheap::Shutdown();

   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      return 1;
   }
}
