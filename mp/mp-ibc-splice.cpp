// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-splice.cpp
//
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

#include "common/environment.h"
#include "common/formatting.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "common/unique.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "mp-algorithms/transfer.h"
#include "tensor/tensor_eigen.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

// The tolerance for the overlap of the left and right boundaries.
double const OverlapTol = 1e-8;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string LeftFilename;
      std::string RightFilename;
      std::string OutputFilename;
      bool Force = false;
      int NLeft = 0;
      int NRight = 0;
      double Tol = 1e-10;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("output,o", prog_opt::value(&OutputFilename), "Output IBC filename [required]")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("tol", prog_opt::value(&Tol), FormatDefault("Tolerance in the squared Frobenius norm of the difference of the center matrix and the translation-invariant fixed point", Tol).c_str())
         ("nleft,n", prog_opt::value(&NLeft), "Number of unit cells to add to the left window [alternative to --tol]")
         ("nright", prog_opt::value(&NRight), "Number of unit cells to add to the right window [if unspecified, will use --nleft if it is specified]")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi-left", prog_opt::value(&LeftFilename), "psi-left")
         ("psi-right", prog_opt::value(&RightFilename), "psi-right")
         ;

      prog_opt::positional_options_description p;
      p.add("psi-left", 1);
      p.add("psi-right", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("psi-right") == 0 || vm.count("output") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi-left> <psi-right> -o <psi-out>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      if (NRight == 0)
         NRight = NLeft;

      if (Verbose > 1)
         std::cout << "Loading wavefunctions..." << std::endl;

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(LeftFilename);
      IBCWavefunction PsiLeft = InPsiLeft->get<IBCWavefunction>();

      pvalue_ptr<MPWavefunction> InPsiRight = pheap::ImportHeap(RightFilename);
      IBCWavefunction PsiRight = InPsiRight->get<IBCWavefunction>();

      CHECK(PsiLeft.right().qshift() == PsiRight.left().qshift());
      QuantumNumber QShift = PsiLeft.right().qshift();

      // Get the middle boundary from PsiLeft in right orthogonal form.
      // (NB LambdaLeft is the Lambda matrix for the translationally invariant
      // fixed point, and will not usually be valid near the window.)
      MatrixOperator LambdaLeft;
      LinearWavefunction BoundaryLeft;
      std::tie(LambdaLeft, BoundaryLeft) = get_right_canonical(PsiLeft.right());
      QuantumNumber LeftQShift = PsiLeft.right_qshift();

      // Get the middle boundary from PsiRight in left orthogonal form.
      MatrixOperator LambdaRight;
      LinearWavefunction BoundaryRight;
      std::tie(BoundaryRight, LambdaRight) = get_left_canonical(PsiRight.left());
      QuantumNumber RightQShift = PsiRight.left_qshift();

      // Get the middle boundary from PsiLeft in LEFT orthogonal form.
      LinearWavefunction BoundaryLeft_LO;
      std::tie(BoundaryLeft_LO, std::ignore) = get_left_canonical(PsiLeft.right());

      if (Verbose > 1)
         std::cout << "Calculating overlap of middle boudaries..." << std::endl;

      std::complex<double> Overlap;
      // Here ULR will be a unitary matrix mapping between the edge basis of BoundaryLeft to BoundaryRight.
      MatrixOperator ULR, RhoLR;
      std::tie(Overlap, ULR, RhoLR) = get_transfer_eigenpair(BoundaryLeft_LO, BoundaryRight, QShift);

      if (std::abs(std::abs(Overlap) - 1.0) > OverlapTol)
      {
         std::cerr << "fatal: The overlap of two middle boundaries is significantly less than 1." << std::endl;
         return 1;
      }

      // Normalize ULR by setting the sum of the singular values of RhoLR to be 1.
      {
         MatrixOperator U, Vh;
         RealDiagonalOperator D;
         SingularValueDecomposition(RhoLR, U, D, Vh);

         RhoLR *= 1.0 / trace(D);
         ULR *= 1.0 / inner_prod(delta_shift(RhoLR, QShift), ULR);
      }

      if (Verbose > 1)
         std::cout << "Building left window..." << std::endl;

      LinearWavefunction WindowLeft;
      MatrixOperator DLeft;
      std::tie(WindowLeft, DLeft) = get_left_canonical(PsiLeft.window());

      if (PsiLeft.window_right_sites() != 0)
      {
         if (Verbose > 1)
            std::cout << "Adding sites to fill the right unit cell of the left window..." << std::endl;

         auto C = BoundaryLeft.begin();
         for (int i = 0; i < PsiLeft.window_right_sites(); ++i)
            ++C;

         MatrixOperator Vh;
         RealDiagonalOperator D;
         while (C != BoundaryLeft.end())
         {
            StateComponent CL = prod(DLeft, delta_shift(*C, LeftQShift));
            std::tie(D, Vh) = OrthogonalizeBasis2(CL);
            WindowLeft.push_back(CL);
            DLeft = D * Vh;
            ++C;
         }
         LeftQShift = delta_shift(LeftQShift, adjoint(QShift));

         // Make DLeft square.
         WindowLeft.set_back(prod(WindowLeft.get_back(), Vh));
         DLeft = adjoint(Vh) * DLeft;
      }
      else
      {
         // Make DLeft square.
         WindowLeft.set_back(prod(WindowLeft.get_back(), PsiLeft.window().RightU()));
         DLeft = adjoint(PsiLeft.window().RightU()) * DLeft;
      }

      if (Verbose > 0)
         std::cout << "NLeft=0 EpsilonLeftSq=" << norm_frob_sq(DLeft - delta_shift(LambdaLeft, LeftQShift)) << std::endl;

      int Iter = 0;
      // If NLeft is unspecified, loop until we reach tolerance, otherwise add NLeft unit cells.
      while (NLeft == 0 ? norm_frob_sq(DLeft - delta_shift(LambdaLeft, LeftQShift)) > Tol : Iter < NLeft)
      {
         ++Iter;

         MatrixOperator Vh;
         RealDiagonalOperator D;
         auto C = BoundaryLeft.begin();
         while (C != BoundaryLeft.end())
         {
            StateComponent CL = prod(DLeft, delta_shift(*C, LeftQShift));
            std::tie(D, Vh) = OrthogonalizeBasis2(CL);
            WindowLeft.push_back(CL);
            DLeft = D * Vh;
            ++C;
         }
         LeftQShift = delta_shift(LeftQShift, adjoint(QShift));

         // Make DLeft square.
         WindowLeft.set_back(prod(WindowLeft.get_back(), Vh));
         DLeft = adjoint(Vh) * DLeft;

         if (Verbose > 0)
            std::cout << "NLeft=" << Iter
                      << " EpsilonLeftSq=" << norm_frob_sq(DLeft - delta_shift(LambdaLeft, LeftQShift)) << std::endl;
      }

      if (norm_frob_sq(DLeft - delta_shift(LambdaLeft, LeftQShift)) > Tol)
         std::cerr << "WARNING: EpsilonLeft=" << norm_frob_sq(DLeft - delta_shift(LambdaLeft, LeftQShift))
                   << " is above tolerance!" << std::endl;

      if (Verbose > 1)
         std::cout << "Building right window..." << std::endl;

      LinearWavefunction WindowRight;
      MatrixOperator DRight;
      std::tie(WindowRight, DRight) = get_left_canonical(PsiRight.window());
      DRight = right_orthogonalize(WindowRight, DRight, Verbose-2);

      if (PsiRight.window_left_sites() != 0)
      {
         if (Verbose > 1)
            std::cout << "Adding sites to fill the left unit cell of the right window..." << std::endl;

         auto C = BoundaryRight.end();
         for (int i = 0; i < PsiRight.window_left_sites(); ++i)
            --C;

         MatrixOperator U;
         RealDiagonalOperator D;
         while (C != BoundaryRight.begin())
         {
            --C;
            StateComponent CR = prod(delta_shift(*C, RightQShift), DRight);
            std::tie(U, D) = OrthogonalizeBasis1(CR);
            WindowRight.push_front(CR);
            DRight = U * D;
         }
         RightQShift = delta_shift(RightQShift, PsiRight.left().qshift());

         // Make DRight square.
         WindowRight.set_front(prod(U, WindowRight.get_front()));
         DRight = DRight * adjoint(U);

      }
      else
      {
         // Make DRight square.
         MatrixOperator U;
         RealDiagonalOperator D;
         StateComponent CR = prod(DRight, WindowRight.get_front());
         std::tie(U, D) = OrthogonalizeBasis1(CR);
         WindowRight.set_front(prod(U, CR));
         DRight = U * D * adjoint(U);
      }

      if (Verbose > 0)
         std::cout << "NRight=0 EpsilonRightSq=" << norm_frob_sq(DRight - delta_shift(LambdaRight, RightQShift)) << std::endl;

      Iter = 0;
      while (NRight == 0 ? norm_frob_sq(DRight - delta_shift(LambdaRight, RightQShift)) > Tol : Iter < NRight)
      {
         ++Iter;

         MatrixOperator U;
         RealDiagonalOperator D;
         auto C = BoundaryRight.end();
         while (C != BoundaryRight.begin())
         {
            --C;
            StateComponent CR = prod(delta_shift(*C, RightQShift), DRight);
            std::tie(U, D) = OrthogonalizeBasis1(CR);
            WindowRight.push_front(CR);
            DRight = U * D;
         }
         RightQShift = delta_shift(RightQShift, PsiRight.left().qshift());

         // Make DRight square.
         WindowRight.set_front(prod(U, WindowRight.get_front()));
         DRight = DRight * adjoint(U);

         if (Verbose > 0)
            std::cout << "NRight=" << Iter
                      << " EpsilonRightSq=" << norm_frob_sq(DRight - delta_shift(LambdaRight, RightQShift)) << std::endl;
      }

      if (norm_frob_sq(DRight - delta_shift(LambdaRight, RightQShift)) > Tol)
         std::cerr << "WARNING: EpsilonRight=" << norm_frob_sq(DRight - delta_shift(LambdaRight, RightQShift))
                   << " is above tolerance!" << std::endl;

      if (Verbose > 1)
         std::cout << "Building spliced window..." << std::endl;

      LinearWavefunction Window;
      for (auto const& C : WindowLeft)
         Window.push_back(delta_shift(C, adjoint(LeftQShift)));

      Window.set_back(prod(Window.get_back(), LambdaLeft * ULR));

      RightQShift = delta_shift(RightQShift, adjoint(QShift));
      for (auto const& C : WindowRight)
         Window.push_back(delta_shift(C, adjoint(RightQShift)));

      MatrixOperator I = MatrixOperator::make_identity(Window.Basis2());
      WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Window), I, Verbose-1);

      if (Verbose > 1)
         std::cout << "Saving wavefunction..." << std::endl;

      IBCWavefunction PsiOut(PsiLeft.left(), PsiWindow, PsiRight.right(),
                             delta_shift(PsiLeft.left_qshift(), adjoint(LeftQShift)),
                             delta_shift(PsiRight.right_qshift(), adjoint(RightQShift)),
                             PsiLeft.window_offset(), PsiLeft.window_left_sites(), PsiRight.window_right_sites());

      // Stream the boundaries, if the input files do.
      PsiOut.set_left_filename(PsiLeft.get_left_filename());
      PsiOut.set_right_filename(PsiRight.get_right_filename());

      MPWavefunction Wavefunction;
      Wavefunction.Wavefunction() = std::move(PsiOut);
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
      Wavefunction.SetDefaultAttributes();

      pvalue_ptr<MPWavefunction> PsiPtr(new MPWavefunction(Wavefunction));
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
