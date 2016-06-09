// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-trispectral.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mp/copyright.h"
#include "common/terminal.h"
#include "common/trace.h"
#include "linearalgebra/eigen.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <vector>
#include <deque>

namespace prog_opt = boost::program_options;
typedef std::complex<double> complex;

complex CalculateSpectral(LinearAlgebra::Vector<double> const& Eigen,
                          LinearAlgebra::Vector<std::complex<double> > const& Components,
                          double Energy,
                          double Frequency,
                          double Broadening)
{
   complex Result = 0;
   for (unsigned i = 0; i < Eigen.size(); ++i)
   {
      //      TRACE(EVec[i])(Eigen[i])(Energy+Frequency-Eigen[i]);
      Result = Result + Components[i] / complex(Energy+Frequency-Eigen[i],Broadening);
   }
   return Result;
}

complex CalculateSpectralBroadeningPropToEnergy(LinearAlgebra::Vector<double> const& Eigen,
                                                LinearAlgebra::Vector<std::complex<double> > const& Components,
                                                double Energy,
                                                double Frequency,
                                                double Broadening)
{
   complex Result = 0;
   for (unsigned i = 0; i < Eigen.size(); ++i)
   {
      //      TRACE(EVec[i])(Eigen[i])(Energy+Frequency-Eigen[i]);
      Result = Result + Components[i] / complex(Energy+Frequency-Eigen[i],
                                                Broadening * (Eigen[i] - Energy + 1e-5));
   }
   return Result;
}

int main(int argc, char** argv)
{
   try
   {
      double Broadening = 0.01;
      double MinFreq;
      double MaxFreq;
      int Count = 21;
      double Interval = 0;
      double GroundstateEnergy = 0;
      bool ShowDos = false;
      bool BroadPropToFreq = false;
      bool BroadPropToEnergy = false;

      std::cout.precision(14);

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("min-freq,m", prog_opt::value<double>(&MinFreq),
          "minimum frequency [required]")
         ("max-freq,x", prog_opt::value<double>(&MaxFreq),
          "maximum frequency [required]")
         ("count,c", prog_opt::value<int>(&Count),
          "Number of points to calculate [default 21]")
         ("interval,i", prog_opt::value<double>(&Interval),
          "alternative to --count, calculate the spectral function at points separated by this interval")
         ("dos,d", prog_opt::bool_switch(&ShowDos),
          "show the normalized density of states, instead of the real and imaginary parts of the Green's function")
         ("Broadening,B", prog_opt::value(&Broadening),
          "width of the Lorentzian broadening of the Green's function [default 0.01]")
         ("GroundstateEnergy,G", prog_opt::value(&GroundstateEnergy),
          "don't load the preamble, instead use this value of the groundstate energy")
         ("broadening-propto-freq", prog_opt::bool_switch(&BroadPropToFreq),
          "Make the broadening proportional to the frequency")
         ("broadening-propto-energy", prog_opt::bool_switch(&BroadPropToEnergy),
          "Make the broadening proportional to the energy")
         ("verbose,v", prog_opt::bool_switch(),
          //("verbose,v", prog_opt::value<std::vector<bool> >(),
          "increase verbosity")
	  ;

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") || vm.count("min-freq") == 0 || vm.count("max-freq") == 0) 
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-trispectral [options] < tridiagonalform\n";
         std::cerr << desc << '\n';
         return 1;
      }

      if (vm.count("count") && vm.count("interval"))
      {
         std::cerr << "fatal: --count and --interval cannot both be specified!\n";
         return 1;
      }

      int const Verbosity = vm["verbose"].as<bool>();
      bool ShowTitles = Verbosity >= 1;

      if (vm.count("GroundstateEnergy") == 0)
         std::cin >> GroundstateEnergy;

      // Load the coefficients of the tridiagonal matrix
      std::deque<double> Alpha, Beta;
      std::vector<std::complex<double> > LanczosLHS;
      int n; double a,b, lvr, lvi;
      while (std::cin >> n >> a >> b >> lvr >> lvi)
      {
         Alpha.push_back(a);
         Beta.push_back(b);
         LanczosLHS.push_back(std::complex<double>(lvr, lvi) * Beta[0]);
      }

      double LanczosNorm = Beta[0];
      Alpha.pop_front();
      Beta.pop_front();

      // M is our matrix
      LinearAlgebra::Matrix<double> M(Alpha.size(), Alpha.size(), 0.0);
      for (unsigned i = 0; i < Alpha.size(); ++i)
      {
         M(i,i) = Alpha[i];
         if (i != Alpha.size()-1)
            M(i+1,i) = M(i,i+1) = Beta[i]; 
      }
      LinearAlgebra::Vector<double> Eigen = DiagonalizeHermitian(M);   // eigenvalues
      LinearAlgebra::Vector<double> EVec = M(LinearAlgebra::all, 0);   // lanczos vector in the eigenbasis

      LinearAlgebra::Vector<std::complex<double> > Components(size1(M));
      // argh, maybe we should add matrix-vector multiply to the matrix lib???
      for (unsigned i = 0; i < size1(M); ++i)
     {
        Components[i] = 0.0;
        for (unsigned j = 0; j < size2(M); ++j)
        {
           Components[i] += M(i,j) * LanczosLHS[j];
        }
        Components[i] = LinearAlgebra::conj(Components[i]) * EVec[i];
     }

      // Set up the std::cout flags
      std::cout.precision(16);
      std::cout.setf(std::ios::showpoint);
      std::cout.setf(std::ios::fixed);

      if (ShowTitles)
      {
         if (ShowDos)
            std::cerr << "Frequency          DOS\n";
         else
            std::cerr << "Frequency          Real           Imag\n";
      }

      double const minus1pi = -1.0 / math_const::pi;
      if (vm.count("interval"))
      {
         // Calculate the spectral function at a fixed interval, rather than a fixed count
         for (double f = MinFreq; f < MaxFreq+Interval/10; f += Interval)
         {
            double ActualBroadening = BroadPropToFreq ? (fabs(Broadening * f)) : Broadening;
            complex x = LanczosNorm*(BroadPropToEnergy ? 
                CalculateSpectralBroadeningPropToEnergy(Eigen, Components, GroundstateEnergy, f, ActualBroadening) 
                                     : CalculateSpectral(Eigen, Components, GroundstateEnergy, f, ActualBroadening));
            std::cout << std::setw(20) << f << "    ";
            if (ShowDos)
               std::cout << std::setw(20) << (minus1pi*x.imag()) << '\n';
            else
               std::cout << std::setw(20) << x.real() << "    "
                         << std::setw(20) << x.imag() << '\n';
         }
      }
      else
      {
         for (int i = 0; i < Count; ++i)
         {
            double f = MinFreq + (MaxFreq-MinFreq) * i / (Count-1);
            double ActualBroadening = BroadPropToFreq ? (fabs(Broadening * f)) : Broadening;
            complex x = LanczosNorm*(BroadPropToEnergy ? 
                CalculateSpectralBroadeningPropToEnergy(Eigen, Components, GroundstateEnergy, f, ActualBroadening) 
                                     : CalculateSpectral(Eigen, Components, GroundstateEnergy, f, ActualBroadening));
            std::cout << std::setw(20) << f << "    ";
            if (ShowDos)
               std::cout << std::setw(20) << (minus1pi*x.imag()) << '\n';
            else
               std::cout << std::setw(20) << x.real() << "    "
                         << std::setw(20) << x.imag() << '\n';
         }
      }

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
