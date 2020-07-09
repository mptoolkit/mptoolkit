// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-expectation.cpp
//
// Copyright (C) 2004-2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

namespace prog_opt = boost::program_options;
using formatting::format_complex;

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
      std::string PsiStr;
      std::string OpStr;
      std::string Psi2Str;
      int Verbose = 0;
      bool Print = false;
      int Coarsegrain = 1;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("real,r", prog_opt::bool_switch(&ShowReal),
          "display only the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImag),
          "display only the imaginary part of the result")
         ("print,p", prog_opt::bool_switch(&Print), "Print the MPO to standard output (use --verbose to see more detail)")
	 ("coarsegrain", prog_opt::value(&Coarsegrain), "coarse-grain N-to-1 sites")
         ("verbose,v", prog_opt_ext::accum_value(&Verbose),
          "Verbose output (use multiple times for more output)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&PsiStr), "psi")
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

      if (vm.count("help") > 0 || vm.count("op") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <operator> [psi2]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // Default formatting a+ib is used unless either --real or --imag is specified
      if (ShowReal || ShowImag)
         ShowDefault = false;

      pvalue_ptr<MPWavefunction> PsiPtr;
      // if we are calculating a mixed expectation value, then we need two wavefunctions so
      // allocate a temporary heap.  Otherwise we can use one heap in read-only mode
      if (vm.count("psi2"))
      {
         mp_pheap::InitializeTempPHeap();
         PsiPtr = pheap::ImportHeap(PsiStr);
      }
      else
      {
         PsiPtr = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);
      }

      UnitCellMPO Op;
      InfiniteLattice Lattice;
      std::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStr);

      CHECK(Op.GetSiteList() == Lattice.GetUnitCell().GetSiteList());

      if (Print)
      {
         print_structure(Op.MPO(), std::cout);
         if (Verbose > 0)
         {
            std::cout << Op.MPO() << '\n';
         }
         if (Verbose > 1)
         {
            SimpleRedOperator x = coarse_grain(Op.MPO());
            std::cout << x << "\n";
         }
         //      std::cout << Op << '\n';
         //std::cout << "\nTransfer matrix:" << construct_transfer_matrix(herm(GenericMPO(Op.MPO())),
         //                                                     GenericMPO(Op.MPO())) << '\n';
      };

      // Check that Op is bosonic, otherwise it is not defined
      CHECK(Op.Commute() == LatticeCommute::Bosonic)("Cannot evaluate non-bosonic operator")(Op.Commute());

      std::complex<double> x; // the expectation value

      if (PsiPtr->is<InfiniteWavefunctionLeft>())
      {
	 if (vm.count("psi2"))
	 {
	    std::cerr << "mp-expectation: fatal: cannot calculate a mixed expectation value of infnite MPS.\n";
	    return 1;
	 }
	 InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

	 // extend Op1 to a multiple of the wavefunction size
	 Op.ExtendToCoverUnitCell(Psi.size() * Coarsegrain);

	 x = expectation(Psi, coarse_grain(Op.MPO(), Coarsegrain));
      }
      else if (PsiPtr->is<FiniteWavefunctionLeft>())
      {
	 FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();
	 FiniteWavefunctionLeft Psi2;
	 if (vm.count("psi2"))
	 {
	    pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(Psi2Str);
	    if (!Psi2Ptr->is<FiniteWavefunctionLeft>())
	    {
	       std::cerr << "mp-expectation: fatal: cannot calculate a mixed expectation value between different types!\n";
	       return 1;
	    }
	    Psi2 = Psi2Ptr->get<FiniteWavefunctionLeft>();
	 }
	 else
	    Psi2 = Psi;

	 Op.ExtendToCoverUnitCell(Psi.size());

	 x = expectation(Psi, Op.MPO(), Psi2);
      }

      if (ShowDefault)
      {
         std::cout << format_complex(x) << '\n';
      }
      else
      {
         if (ShowReal)
            std::cout << x.real() << "   ";
         if (ShowImag)
            std::cout << x.imag();
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
