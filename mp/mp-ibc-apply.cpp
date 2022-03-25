// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-apply.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
#include "wavefunction/mpwavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "lattice/infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "lattice/unitcell-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OpStrLeft;
      std::string OpStrWindow;
      std::string InputFile;
      std::string OutputFile;
      bool Force = false;
      bool Normalize = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("left,l", prog_opt::value(&OpStrLeft),
          "operator for the left boundary")
         ("window,w", prog_opt::value(&OpStrWindow),
          "operator for the window")
         ("force,f", prog_opt::bool_switch(&Force),
          "allow overwriting the output file, if it already exists")
         ("normalize", prog_opt::bool_switch(&Normalize),
          "normalize the output wavefunction")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
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

      if (vm.count("help") > 0 || vm.count("psi2") < 1 || (vm.count("left") < 1 && vm.count("window") < 1))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <input-psi> <output-psi>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFile == OutputFile)
         PsiPtr = pheap::OpenPersistent(InputFile.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFile, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFile);
      }

      IBCWavefunction Psi = PsiPtr->get<IBCWavefunction>();
      InfiniteWavefunctionLeft PsiLeft = Psi.Left;
      InfiniteWavefunctionRight PsiRight = Psi.Right;
      MatrixOperator Vh = MatrixOperator::make_identity(Psi.Left.Basis1());

      if (vm.count("left") > 0) {
         if (Verbose > 0)
            std::cout << "Applying operator to left semi-infinite boundary..." << std::endl;

         InfiniteLattice Lattice;
         ProductMPO StringOp;
         std::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(OpStrLeft);

         if (PsiLeft.size() != StringOp.size())
         {
            int Size = statistics::lcm(PsiLeft.size(), StringOp.size());
            if (PsiLeft.size() < Size)
            {
               std::cerr << "Warning: extending left wavefunction size to lowest common multiple, which is "
                         << Size << " sites." << std::endl;
            }
            PsiLeft = repeat(PsiLeft, Size / PsiLeft.size());
            StringOp = repeat(StringOp, Size / StringOp.size());
         }

         LinearWavefunction PsiLeftLinear = get_left_canonical(PsiLeft).first;

         if (Verbose > 0)
            std::cout << "Applying operator..." << std::endl;

         // apply the string operator
         ProductMPO::const_iterator MI = StringOp.begin();
         for (LinearWavefunction::iterator I = PsiLeftLinear.begin(); I != PsiLeftLinear.end(); ++I, ++MI)
         {
            (*I) = aux_tensor_prod(*MI, *I);
         }

         // TODO: At the moment, the function which orthgonalizes the state
         // causes issues, so we assume that the operator does not effect the
         // orthogonality of the unit cell and use ConstructFromOrthogonal.
#if 0
         PsiLeft = InfiniteWavefunctionLeft::Construct(PsiLeftLinear, PsiLeft.qshift(), Verbose);
#else
         PsiLeft = InfiniteWavefunctionLeft::ConstructFromOrthogonal(PsiLeftLinear, PsiLeft.lambda_r(), PsiLeft.qshift(), Vh, Verbose);
#endif
         Vh = delta_shift(Vh, adjoint(PsiLeft.qshift()));

         if (Verbose > 0)
            std::cout << "Finished applying operator to left semi-infinite boundary..." << std::endl;
      }

      // The number of sites to incorporate from the left/right boundaries, if necessary.
      int SitesLeft = 0;
      int SitesRight = 0;

      int NewOffset = Psi.window_offset();

      WavefunctionSectionLeft PsiWindow = Psi.Window;

      LinearWavefunction PsiWindowLinear;
      MatrixOperator Lambda;

      // For an empty window, we cannot use get_left_canonical.
      if (Psi.window_size() > 0)
      {
         std::tie(PsiWindowLinear, Lambda) = get_left_canonical(PsiWindow);
         PsiWindowLinear.set_back(prod(PsiWindowLinear.get_back(), Lambda));
         PsiWindowLinear.set_front(prod(Vh, PsiWindowLinear.get_front()));
      }
      else
      {
         PsiWindowLinear = LinearWavefunction();
         Lambda = Vh * PsiWindow.LeftU() * PsiWindow.lambda_r() * PsiWindow.RightU();
      }

      if (vm.count("window") > 0) {
         if (Verbose > 0)
            std::cout << "Applying operator to window..." << std::endl;

         InfiniteLattice Lattice;
         UnitCellMPO Op;
         std::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStrWindow);

         // Calculate the number of sites that we need to incorporate to the
         // window from the left and right semi-infinite boundaries.
         SitesLeft = std::max(Psi.window_offset() - Op.offset(), 0);
         SitesRight = std::max(Op.size()+Op.offset() - Psi.window_size()-Psi.window_offset(), 0);

         // Ensure that the window starts and ends at the operator unit cell boundaries.
         SitesLeft += (Psi.Left.size() - Psi.WindowLeftSites) % Op.unit_cell_size();
         SitesRight += (Psi.Right.size() - Psi.WindowRightSites) % Op.unit_cell_size();

         // The new window offset.
         NewOffset = Psi.window_offset() - SitesLeft;

         // Incorporate extra sites from the left boundary.
         if (Verbose > 0)
            std::cout << "Incorporating " << SitesLeft
                      << " sites from the left boundary..." << std::endl;

         auto CLeft = PsiLeft.end();

         for (int i = 0; i < Psi.WindowLeftSites; ++i)
            --CLeft;

         for (int i = 0; i < SitesLeft; ++i)
         {
            if (CLeft == PsiLeft.begin())
            {
               inplace_qshift(PsiLeft, PsiLeft.qshift());
               CLeft = PsiLeft.end();
            }
            --CLeft;
            PsiWindowLinear.push_front(*CLeft);
         }

         // If the initial window had no sites and we just added sites from the
         // left, incorporate the lambda matrix now; otherwise, we incorporate it
         // below after adding from the right.
         if (SitesLeft > 0 && Psi.window_size() == 0)
            PsiWindowLinear.set_back(prod(PsiWindowLinear.get_back(), Lambda));

         // Incorporate extra sites from the right boundary.
         if (Verbose > 0)
            std::cout << "Incorporating " << SitesRight
                      << " sites from the right boundary..." << std::endl;

         auto CRight = PsiRight.begin();

         for (int i = 0; i < Psi.WindowRightSites; ++i)
            ++CRight;

         for (int i = 0; i < SitesRight; ++i)
         {
            if (CRight == PsiRight.end())
            {
               inplace_qshift(PsiRight, adjoint(PsiRight.qshift()));
               CRight = PsiRight.begin();
            }
            PsiWindowLinear.push_back(*CRight);
            ++CRight;
         }

         if (SitesLeft == 0 && Psi.window_size() == 0)
            PsiWindowLinear.set_front(prod(Lambda, PsiWindowLinear.get_front()));

         Op.ExtendToCover(PsiWindowLinear.size(), NewOffset);

         if (Verbose > 0)
            std::cout << "Applying operator..." << std::endl;

         BasicFiniteMPO::const_iterator MI = Op.MPO().begin();

         for (auto I = PsiWindowLinear.begin(); I != PsiWindowLinear.end(); ++I, ++MI)
            (*I) = aux_tensor_prod(*MI, *I);

         // Handle the case where the MPO has a nontrivial quantum number shift.
         // (This only works when Op.Basis1.size() == 1).
         inplace_qshift(PsiLeft, Op.qn1());

         if (Verbose > 0)
            std::cout << "Finished applying operator to window..." << std::endl;
      }

      if (Psi.window_size() + SitesLeft + SitesRight > 0)
      {
         MatrixOperator Identity = MatrixOperator::make_identity(PsiWindowLinear.Basis2());
         PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), Identity, Verbose);

         if (Normalize)
         {
            if (Verbose > 0)
               std::cout << "Normalizing wavefunction..." << std::endl;

            std::tie(PsiWindowLinear, Lambda) = get_left_canonical(PsiWindow);
            Lambda *= 1.0 / norm_frob(Lambda);
            PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), Lambda, Verbose);
         }
      }
      else
      {
         if (Normalize)
         {
            if (Verbose > 0)
               std::cout << "Normalizing wavefunction..." << std::endl;
            Lambda *= 1.0 / norm_frob(Lambda);
         }
         PsiWindow = WavefunctionSectionLeft(Lambda);
      }

      IBCWavefunction PsiNew;
      PsiNew = IBCWavefunction(PsiLeft, PsiWindow, PsiRight, NewOffset,
                               (Psi.WindowLeftSites + SitesLeft) % PsiLeft.size(),
                               (Psi.WindowRightSites + SitesRight) % PsiRight.size());

      PsiPtr.mutate()->Wavefunction() = PsiNew;

      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
