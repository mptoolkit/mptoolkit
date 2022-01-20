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
#include "lattice/unitcell-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;


int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string OpStr;
      std::string InputFile;
      std::string OutputFile;
      bool Force = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("force,f", prog_opt::bool_switch(&Force),
          "allow overwriting the output file, if it already exists")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("op", prog_opt::value(&OpStr), "op")
         ("psi1", prog_opt::value(&InputFile), "psi1")
         ("psi2", prog_opt::value(&OutputFile), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("op", 1);
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
         std::cerr << "usage: " << basename(argv[0]) << " [options] <operator> <input-psi> <output-psi>\n";
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

      InfiniteLattice Lattice;
      UnitCellMPO Op;
      std::tie(Op, Lattice) = ParseUnitCellOperatorAndLattice(OpStr);

      LinearWavefunction PsiWindow;
      MatrixOperator Lambda;

      // For an empty window, we cannot use get_left_canonical.
      if (Psi.window_size() > 0)
      {
         std::tie(PsiWindow, Lambda) = get_left_canonical(Psi.Window);
         PsiWindow.set_back(prod(PsiWindow.get_back(), Lambda));
      }
      else
      {
         PsiWindow = LinearWavefunction();
         Lambda = Psi.Window.LeftU() * Psi.Window.lambda_r() * Psi.Window.RightU();
      }

      InfiniteWavefunctionLeft PsiLeft = Psi.Left;
      InfiniteWavefunctionRight PsiRight = Psi.Right;

      // Calculate the number of sites that we need to incorporate to the
      // window from the left and right semi-infinite boundaries.
      int SitesLeft = std::max(Psi.window_offset() - Op.offset(), 0);
      int SitesRight = std::max(Op.size()+Op.offset() - Psi.window_size()-Psi.window_offset(), 0);

      // Ensure that the window starts and ends at the operator unit cell boundaries.
      SitesLeft += (Psi.Left.size() - Psi.WindowLeftSites) % Op.unit_cell_size();
      SitesRight += (Psi.Right.size() - Psi.WindowRightSites) % Op.unit_cell_size();

      // The new window offset.
      int NewOffset = Psi.window_offset() - SitesLeft;

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
         PsiWindow.push_front(*CLeft);
      }

      // If the initial window had no sites and we just added sites from the
      // left, incorporate the lambda matrix now; otherwise, we incorporate it
      // below after adding from the right.
      if (SitesLeft > 0 && Psi.window_size() == 0)
         PsiWindow.set_back(prod(PsiWindow.get_back(), Lambda));

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
         PsiWindow.push_back(*CRight);
         ++CRight;
      }

      if (SitesLeft == 0 && Psi.window_size() == 0)
         PsiWindow.set_front(prod(Lambda, PsiWindow.get_front()));

      Op.ExtendToCover(PsiWindow.size(), NewOffset);

      if (Verbose > 0)
         std::cout << "Applying operator..." << std::endl;

      BasicFiniteMPO::const_iterator MI = Op.MPO().begin();

      for (auto I = PsiWindow.begin(); I != PsiWindow.end(); ++I, ++MI)
         (*I) = aux_tensor_prod(*MI, *I);

      MatrixOperator Identity = MatrixOperator::make_identity(PsiWindow.Basis2());
      WavefunctionSectionLeft PsiWindowCanonical = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindow), Identity, Verbose);

      // Handle the case where the MPO has a nontrivial quantum number shift.
      // (This only works when Op.Basis1.size() > 1).
      inplace_qshift(PsiLeft, Op.qn1());

      IBCWavefunction PsiNew;
      PsiNew = IBCWavefunction(PsiLeft, PsiWindowCanonical, PsiRight, NewOffset,
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
