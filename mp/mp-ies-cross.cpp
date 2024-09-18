// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-ies-cross.cpp
//
// Copyright (C) 2015-2024 Ian McCulloch <ian@qusim.net>
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

#include "wavefunction/mpwavefunction.h"
#include "mp-algorithms/transfer.h"
#include "lattice/latticesite.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>
#include "common/formatting.h"
#include <tuple>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace prog_opt = boost::program_options;
using formatting::format_complex;

// returns true if Name exists and is a regular file
bool FileExists(std::string const& Name)
{
   struct stat buf;
   return stat(Name.c_str(), &buf) != -1 && S_ISREG(buf.st_mode);
}

// The eigenvalues of the operator, quantum number sector, magnitude and
// phase angle
struct EValueRecord
{
   QuantumNumber q;
   std::complex<double> Eigenvalue;
};

bool CompareMagnitude(EValueRecord const& x, EValueRecord const& y)
{
   return std::abs(x.Eigenvalue) > std::abs(y.Eigenvalue);
}

bool CompareQMagnitude(EValueRecord const& x, EValueRecord const& y)
{
   return x.q < y.q || (x.q == y.q && std::abs(x.Eigenvalue) > std::abs(y.Eigenvalue));
}

void PrintHeading(std::ostream& out, bool ShowRadians)
{
   out << "#Number #Sector     #Degen  #Angle" << (ShowRadians ? "_rad " : "_deg ")
       <<  "          #Eigenvalue             #Energy\n";
}

void PrintFormat(std::ostream& out, int n, EValueRecord const& x, bool ShowRadians)
{
   std::string SectorStr = boost::lexical_cast<std::string>(x.q);
   out << std::left << std::setw(7) << n << ' ';
   out << std::left << std::setw(11) << SectorStr << ' ';
   out << std::left << std::setw(6) << degree(x.q) << ' ';
   double Arg = std::arg(x.Eigenvalue);
   if (!ShowRadians)
      Arg *= 180.0 / math_const::pi;
   out << std::setw(20) << std::right << std::fixed << Arg << "  ";
   out << std::setw(20) << std::scientific << std::abs(x.Eigenvalue) << ' ';
   out << std::setw(20) << std::fixed << (-log(std::abs(x.Eigenvalue)));
   out << std::setfill(' ') << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      std::string Psi1Str;
      std::string Psi2Str;
      double Tol = 1E-14;
      int Verbose = 0;
      bool ShowRadians = false;
      bool Quiet = false;
      bool SplitBySectors = false;
      std::string QSector;
      std::string SplitOutputPrefix;
      bool Entropy = false;
      std::vector<double> Renyi;
      bool Conj = false;
      bool Reflect = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show column headings")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "display the angle in radians instead of degrees (also for the --gauge angle, if specified)")
         ("qsector,q", prog_opt::value(&QSector), "Quantum number sector of the transfer operator [default scalar sector]")
         ("entropy,e", prog_opt::bool_switch(&Entropy), "Calculate the von Neumann entropy")
         ("reyni,a", prog_opt::value(&Renyi), "Calculate the Renyi entropy for this value of alpha (can be used more than once)")
         ("conj", prog_opt::bool_switch(&Conj), "complex conjugate psi1")
         ("reflect", prog_opt::bool_switch(&Reflect), "reflect psi1")
         ("sectors", prog_opt::bool_switch(&SplitBySectors), "sort output by quantum number sector")
         ("sectors-prefix", prog_opt::value(&SplitOutputPrefix),
          "Write the output to a separate file for each quantum number sector prefix.sector.dat")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value<std::string>(&Psi1Str), "psi1")
         ("psi2", prog_opt::value<std::string>(&Psi2Str), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi1", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi1> <psi2>\n";
         std::cerr << desc << '\n';
         std::cerr << "This tool constructs the entanglement spectrum of  |Ψ₁><Ψ₂|.\n"
                   << "The phase angles are only determined up to a global phase (which\n"
                   << "arises from the global phase of the eigenmatrix of the symmetry\n"
                   << "operator).\n"
                   << "\nBy default the entanglement spectrum is printed to standard output,\n"
                   << "ordered by eigenvalue.  Alternatively, the spectrum can be displayed\n"
                   << "for each quantum number sector separately, with the --sectors option.\n"
                   << "The --sectors-prefix=<fileprefix> option implies --sectors, and writes\n"
                   << "each sector to a file of the form fileprefix-q<sector>.dat\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (vm.count("sectors-prefix"))
         SplitBySectors = true;

      // Load the wavefunction
      pvalue_ptr<MPWavefunction> Psi1Ptr = pheap::OpenPersistent(Psi1Str, mp_pheap::CacheSize(), true);
      pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(Psi2Str);

      InfiniteWavefunctionLeft Psi1 = Psi1Ptr->get<InfiniteWavefunctionLeft>();
      InfiniteWavefunctionLeft Psi2 = Psi2Ptr->get<InfiniteWavefunctionLeft>();
      auto q = QuantumNumber(Psi1.GetSymmetryList());
      if (vm.count("qsector"))
         q = QuantumNumber(Psi1.GetSymmetryList(), QSector);

      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi1..." << std::endl;
         inplace_reflect(Psi1);
      }
      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi1..." << std::endl;
         inplace_conj(Psi1);
      }

      // Get the matrix
      std::complex<double> e;
      MatrixOperator L, R;
      std::tie(e, L, R) = get_transfer_eigenpair(Psi1, Psi2, q);

      MatrixOperator v = L*R;
      v *= 1.0 / trace(v);

      // Now we want the eigenvalues of v
      // These will be of the form a * exp(i*theta) where a is the density matrix eigenvalue and
      // theta is the phase (up to an overall phase ambiguity, since it is a constant)

      std::vector<EValueRecord> EValues;
      std::vector<std::complex<double>> RenyiEntropy(Renyi.size(), 0.0);
      std::complex<double> vnEntropy = 0.0;

      for (std::size_t i = 0; i < v.Basis1().size(); ++i)
      {
         if (iterate_at(v.data(), i, i))
         {
            LinearAlgebra::Vector<std::complex<double>> Eig =
               LinearAlgebra::EigenvaluesComplex(v(i,i));

            for (unsigned j = 0; j < size(Eig); ++j)
            {
               EValues.push_back({v.Basis1()[i], Eig[j]});
               vnEntropy += -Eig[j] * std::log(Eig[j]);
               for (int k = 0; k < Renyi.size(); ++k)
               {
                  RenyiEntropy[k] += std::pow(Eig[j], Renyi[k]);
               }
            }
         }
      }

      // sort
      if (SplitBySectors)
         std::sort(EValues.begin(), EValues.end(), &CompareQMagnitude);
      else
         std::sort(EValues.begin(), EValues.end(), &CompareMagnitude);

      // output

      if (SplitBySectors)
      {
         QuantumNumber Sector;
         int k = 0;
         int kOffset = 0;
         std::ofstream OutFile;
         std::ostream& Out = SplitOutputPrefix.empty() ? std::cout : OutFile;
         if (SplitOutputPrefix.empty() && !Quiet)
         {
            print_preamble(Out, argc, argv);
            Out << "#Mixed transfer matrix eigenvalue is " << format_complex(e) << '\n';
         }

         while (k < int(EValues.size()))
         {
            if (Sector.is_null() || EValues[k].q != Sector)
            {
               kOffset = k;
               Sector = EValues[k].q;
               if (!SplitOutputPrefix.empty())
               {
                  OutFile.close();
                  OutFile.open(SplitOutputPrefix + "-q" + boost::lexical_cast<std::string>(EValues[k].q) + ".dat",
                               std::ios::out | std::ios::trunc);
                  OutFile.precision(getenv_or_default("MP_PRECISION", 14));
                  if (!Quiet)
                  {
                     print_preamble(Out, argc, argv);
                     Out << "#Eigenvalue of operator is " << format_complex(e) << '\n';
                     Out << "#Sector: " << Sector << '\n';
                     PrintHeading(Out, ShowRadians);
                  }
               }
               else
               {
                  if (!Quiet)
                  {
                     Out << "\n#Sector: " << Sector << '\n';
                     PrintHeading(Out, ShowRadians);
                  }
               }
            }
            PrintFormat(Out, k-kOffset, EValues[k], ShowRadians);
            ++k;
         }
      }
      else
      {
         if (!Quiet)
         {
            print_preamble(std::cout, argc, argv);
            std::cout << "#Eigenvalue of operator is " << format_complex(e) << '\n';
            PrintHeading(std::cout, ShowRadians);
         }
         int k = 0; // eigenvalue number
         for (auto const& e : EValues)
         {
            PrintFormat(std::cout, k++, e, ShowRadians);
         }
      }

      if (Entropy || !Renyi.empty())
         std::cout << '\n';

      if (Entropy)
         std::cout << "Entropy = " << format_complex(vnEntropy) << '\n';

      if (!Renyi.empty())
      {
         for (int i = 0; i < Renyi.size(); ++i)
         {
            std::cout << "Renyi entropy (alpha=" << Renyi[i] << ") = " << format_complex((1.0 / (1.0 - Renyi[i])) * std::log(RenyiEntropy[i])) << '\n';
         }
      }

      pheap::Shutdown();
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
