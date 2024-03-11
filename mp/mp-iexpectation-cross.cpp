// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-iexpectation-cross.cpp
//
// Copyright (C) 2012-2021 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "common/prog_options.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell-parser.h"
#include "wavefunction/operator_actions.h"
#include "mp-algorithms/transfer.h"

namespace prog_opt = boost::program_options;
using formatting::format_complex;

// Lambda si the transfer matrix eigenvalue (per wavefunction unit cell size),
// e is the expectation value
struct expectation_result
{
   std::complex<double> Lambda;
   std::complex<double> e;
};

expectation_result
expectation_cross(InfiniteWavefunctionLeft const& Psi1, BasicFiniteMPO const& Op, InfiniteWavefunctionLeft const& Psi2,
                  QuantumNumber const& q)
{
   std::complex<double> e;
   MatrixOperator Left, Right;
   std::tie(e, Left, Right) = get_transfer_eigenpair(Psi1, Psi2, q);

   MatrixOperator X = Left;
   X = inject_left(X, Psi1, Op, Psi2);

   auto c = inner_prod(delta_shift(Right, Psi1.qshift()), X);

   CHECK_EQUAL(Op.size() % Psi1.size(), 0);
   c *= std::pow(e, -Op.size() / Psi1.size());
   return {e, c};
}

void DisplayHeading(bool ShowReal, bool ShowImag)
{
   // Output the heading
   std::cout << "#i #j";
   if (ShowReal)
      std::cout << " #real";
   if (ShowImag)
      std::cout << " #imag";
   std::cout << '\n';
}

void
Display(std::complex<double> x, int s1, int s2, bool ShowReal, bool ShowImag)
{
   std::cout << s1 << "    " << s2 << "   ";
   if (ShowReal)
      std::cout << x.real() << "   ";
   if (ShowImag)
      std::cout << x.imag();
   std::cout << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      bool ShowReal = false, ShowImag = false;
      bool ShowDefault = true;
      std::string Psi1Str;
      std::string OpStr;
      std::string Psi2Str;
      std::string CacheDirectory;
      std::string Sector;
      bool NoCache = false;
      bool Quiet = false;
      int Verbose = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("real,r", prog_opt::bool_switch(&ShowReal), "display only the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImag), "display only the imaginary part of the result")
         ("quantumnumber,q", prog_opt::value(&Sector), "quantum number sector of the transfer matrix")
         ("nocache", prog_opt::bool_switch(&NoCache), "don't cache the transfer matrix eigenvectors")
         ("cachedir", prog_opt::value(&CacheDirectory), "use this cache directory")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose), "Verbose output (use multiple times for more output)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&Psi1Str), "psi")
         ("op", prog_opt::value(&OpStr), "op")
         ("psi2", prog_opt::value(&Psi2Str), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("op", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("psi2") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <operator> <psi2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // Default formatting a+ib is used unless either --real or --imag is specified
      if (ShowReal || ShowImag)
         ShowDefault = false;

      pvalue_ptr<MPWavefunction> Psi1Ptr;
      pvalue_ptr<MPWavefunction> Psi2Ptr;
      mp_pheap::InitializeTempPHeap();
      Psi1Ptr = pheap::ImportHeap(Psi1Str);
      Psi2Ptr = pheap::ImportHeap(Psi2Str);

      UnitCellMPO Op;
      InfiniteLattice Lattice;
      std::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStr);

      CHECK(Op.GetSiteList() == Lattice.GetUnitCell().GetSiteList());

      // Check that Op is bosonic, otherwise it is not defined
      CHECK(Op.Commute() == LatticeCommute::Bosonic)("Cannot evaluate non-bosonic operator")(Op.Commute());

      if (!Psi1Ptr->is<InfiniteWavefunctionLeft>())
      {
         std::cerr << basename(argv[0]) << ": fatal: expected an infinite wavefunction for psi1.\n";
         return 1;
      }
      if (!Psi2Ptr->is<InfiniteWavefunctionLeft>())
      {
         std::cerr << basename(argv[0]) << ": fatal: expected an infinite wavefunction for psi2.\n";
         return 1;
      }

      InfiniteWavefunctionLeft Psi1 = Psi1Ptr->get<InfiniteWavefunctionLeft>();
      InfiniteWavefunctionLeft Psi2 = Psi2Ptr->get<InfiniteWavefunctionLeft>();

      // extend Op1 to a multiple of the wavefunction size
      Op.ExtendToCoverUnitCell(Psi1.size());

      auto q = QuantumNumber(Psi1.GetSymmetryList(), Sector);

      expectation_result R = expectation_cross(Psi1, Op.MPO(), Psi2, q);

      if (Verbose > 0)
      {
         std::cerr << "Transfer matrix eigenvalue is " << format_complex(R.Lambda) << " per unit cell size " << Psi1.size() << '\n';
      }

      if (ShowDefault)
      {
         std::cout << format_complex(R.e) << '\n';
      }
      else
      {
         if (ShowReal)
            std::cout << R.e.real() << "   ";
         if (ShowImag)
            std::cout << R.e.imag();
         std::cout << '\n';
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
