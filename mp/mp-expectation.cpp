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

// TODO:
// Move this function to the correct place.
// Handle the case where the window is offset by a number of sites not
// divisible by the unit cell size.
std::complex<double>
expectation(IBCWavefunction const& Psi,
            UnitCellMPO Op)
{
   LinearWavefunction PsiLinear;
   MatrixOperator Lambda;

   // For an empty window, we cannot use get_left_canonical.
   if (Psi.window_size() > 0)
   {
      std::tie(PsiLinear, Lambda) = get_left_canonical(Psi.Window);
      PsiLinear.set_back(prod(PsiLinear.get_back(), Lambda));
   }
   else
   {
      PsiLinear = LinearWavefunction();
      Lambda = Psi.Window.lambda_r();
   }

   // Add extra sites to the left, if needed.
   LinearWavefunction PsiLeft;
   RealDiagonalOperator LambdaLeft;
   std::tie(PsiLeft, LambdaLeft) = get_left_canonical(Psi.Left);

   int SitesLeft = std::max(Psi.window_offset() - Op.offset(), 0);

   LinearWavefunction::const_iterator CLeft = PsiLeft.end();
   for (int i = 0; i < SitesLeft; ++i)
   {
      --CLeft;
      PsiLinear.push_front(*CLeft);
      if (CLeft == PsiLeft.begin())
         CLeft = PsiLeft.end();
   }

   if (Psi.window_size() == 0 && SitesLeft > 0)
      PsiLinear.set_back(prod(PsiLinear.get_back(), Lambda));

   // Add extra sites to the right, if needed.
   LinearWavefunction PsiRight;
   RealDiagonalOperator LambdaRight;
   std::tie(LambdaRight, PsiRight) = get_right_canonical(Psi.Right);

   int SitesRight = std::max(Op.size()+Op.offset() - Psi.window_size()-Psi.window_offset(), 0);

   LinearWavefunction::const_iterator CRight = PsiRight.begin();
   for (int i = 0; i < SitesRight; ++i)
   {
      PsiLinear.push_back(*CRight);
      ++CRight;
      if (CRight == PsiRight.end())
         CRight = PsiRight.begin();
   }

   if (Psi.window_size() == 0 && SitesLeft == 0)
      PsiLinear.set_front(prod(Lambda, PsiLinear.get_front()));

   Op.ExtendToCover(PsiLinear.size(), Psi.window_offset()-SitesLeft);

   BasicFiniteMPO M = Op.MPO();

   MatrixOperator I = MatrixOperator::make_identity(PsiLinear.Basis1());
   StateComponent E(M.Basis1(), I.Basis1(), I.Basis2());
   E[0] = I;

   LinearWavefunction::const_iterator C = PsiLinear.begin();
   BasicFiniteMPO::const_iterator W = M.begin();
   while (C != PsiLinear.end())
   {
      E = contract_from_left(*W, herm(*C), E, *C);
      ++C, ++W;
   }

   return trace(E[0]);
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
            std::cerr << "mp-expectation: fatal: cannot calculate a mixed expectation value of infinite MPS.\n"
            "Use mp-iexpectation-cross instead.";
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
      else if (PsiPtr->is<IBCWavefunction>())
      {
         if (vm.count("psi2"))
         {
            std::cerr << "mp-expectation: fatal: cannot calculate a mixed expectation value of IBC wavefunctions.\n";
            return 1;
         }

         IBCWavefunction Psi = PsiPtr->get<IBCWavefunction>();

         x = expectation(Psi, Op);
      }
      else
      {
         std::cerr << "mp-expectation: fatal: unknown wavefunction type.\n";
         return 1;
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
