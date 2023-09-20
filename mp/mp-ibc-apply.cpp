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
      std::string InputFilename;
      std::string OutputFilename;
      bool Force = false;
      bool Normalize = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("left,l", prog_opt::value(&OpStrLeft), "Operator for the left boundary")
         ("window,w", prog_opt::value(&OpStrWindow), "Operator for the window")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("normalize", prog_opt::bool_switch(&Normalize), "Normalize the output wavefunction")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi1", prog_opt::value(&InputFilename), "psi1")
         ("psi2", prog_opt::value(&OutputFilename), "psi2")
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
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi-in> <psi-out>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      pvalue_ptr<MPWavefunction> PsiPtr;
      if (InputFilename == OutputFilename)
         PsiPtr = pheap::OpenPersistent(InputFilename.c_str(), mp_pheap::CacheSize());
      else
      {
         pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);
         PsiPtr = pheap::ImportHeap(InputFilename);
      }

      IBCWavefunction Psi = PsiPtr->get<IBCWavefunction>();
      InfiniteWavefunctionLeft PsiLeft = Psi.left();
      QuantumNumbers::QuantumNumber LeftQShift = Psi.left_qshift();
      InfiniteWavefunctionRight PsiRight = Psi.right();
      QuantumNumbers::QuantumNumber RightQShift = Psi.right_qshift();

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
         PsiLeft = InfiniteWavefunctionLeft::ConstructFromOrthogonal(PsiLeftLinear, PsiLeft.qshift(), PsiLeft.lambda_r(), Verbose);
#endif

         if (Verbose > 0)
            std::cout << "Finished applying operator to left semi-infinite boundary..." << std::endl;
      }

      // The number of sites to incorporate from the left/right boundaries, if necessary.
      int SitesLeft = 0;
      int SitesRight = 0;

      int NewOffset = Psi.window_offset();

      WavefunctionSectionLeft PsiWindow = Psi.window();

      LinearWavefunction PsiWindowLinear;
      MatrixOperator Lambda;

      // For an empty window, we cannot use get_left_canonical.
      if (Psi.window_size() > 0)
      {
         std::tie(PsiWindowLinear, Lambda) = get_left_canonical(PsiWindow);
         PsiWindowLinear.set_back(prod(PsiWindowLinear.get_back(), Lambda));
      }
      else
      {
         PsiWindowLinear = LinearWavefunction();
         Lambda = PsiWindow.LeftU() * PsiWindow.lambda_r() * PsiWindow.RightU();
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
         SitesLeft += (PsiLeft.size() - Psi.window_left_sites()) % Op.unit_cell_size();
         SitesRight += (PsiRight.size() - Psi.window_right_sites()) % Op.unit_cell_size();

         // The new window offset.
         NewOffset = Psi.window_offset() - SitesLeft;

         // Incorporate extra sites from the left boundary.
         if (Verbose > 0)
            std::cout << "Incorporating " << SitesLeft
                      << " sites from the left boundary..." << std::endl;

         auto CLeft = PsiLeft.end();

         for (int i = 0; i < Psi.window_left_sites(); ++i)
            --CLeft;

         for (int i = 0; i < SitesLeft; ++i)
         {
            if (CLeft == PsiLeft.begin())
            {
               LeftQShift = delta_shift(LeftQShift, PsiLeft.qshift());
               CLeft = PsiLeft.end();
            }
            --CLeft;
            PsiWindowLinear.push_front(delta_shift(*CLeft, LeftQShift));
         }

         if (CLeft == PsiLeft.begin())
            LeftQShift = delta_shift(LeftQShift, PsiLeft.qshift());

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

         for (int i = 0; i < Psi.window_right_sites(); ++i)
            ++CRight;

         for (int i = 0; i < SitesRight; ++i)
         {
            if (CRight == PsiRight.end())
            {
               RightQShift = delta_shift(RightQShift, adjoint(PsiRight.qshift()));
               CRight = PsiRight.begin();
            }
            PsiWindowLinear.push_back(delta_shift(*CRight, RightQShift));
            ++CRight;
         }

         if (CRight == PsiRight.end())
            RightQShift = delta_shift(RightQShift, adjoint(PsiRight.qshift()));

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
         LeftQShift = delta_shift(LeftQShift, Op.qn1());

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
      PsiNew = IBCWavefunction(PsiLeft, PsiWindow, PsiRight, LeftQShift, RightQShift, NewOffset,
                               (Psi.window_left_sites() + SitesLeft) % PsiLeft.size(),
                               (Psi.window_right_sites() + SitesRight) % PsiRight.size());

      // Stream the boundaries, if the input file does.
      // UNLESS we modify the left boundary.
      if (!Psi.get_left_filename().empty() && vm.count("left") == 0)
         PsiNew.set_left_filename(Psi.get_left_filename());

      PsiNew.set_right_filename(Psi.get_right_filename());

      PsiPtr.mutate()->Wavefunction() = PsiNew;

      PsiPtr.mutate()->AppendHistoryCommand(EscapeCommandline(argc, argv));
      PsiPtr.mutate()->SetDefaultAttributes();

      if (Verbose > 0)
         std::cout << "Finished." << std::endl;

      pheap::ShutdownPersistent(PsiPtr);
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      if (e.Why == "File exists")
         std::cerr << "Note: Use --force (-f) option to overwrite." << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      pheap::Cleanup();
      return 1;
   }
}
