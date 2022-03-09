// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-overlap.cpp
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

// The tolerance for the left/right boundary overlaps for the overlap of two general IBCs.
double const OverlapTol = 1e-12;

// The tolerance of the trace of the left/right boundary eigenvectors for
// fixing the phase when calculating the overlap of two general IBCs.
double const TraceTol = 1e-8;

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

class ConstIBCIterator
{
   public:
      ConstIBCIterator(IBCWavefunction const& Psi_, int Index_)
         : Psi(Psi_), Index(Index_)
      {
         PsiLeft = Psi.Left;
         PsiRight = Psi.Right;

         WindowLeftIndex = Psi.window_offset();
         WindowRightIndex = Psi.window_size() + Psi.window_offset() - 1;

         if (Index < WindowLeftIndex)
         {
            int IndexDiff = WindowLeftIndex - Index;

            for (int i = 0; i < (Psi.WindowLeftSites + IndexDiff) / PsiLeft.size(); ++i)
               inplace_qshift(PsiLeft, PsiLeft.qshift());

            C = PsiLeft.end();
            for (int i = 0; i < (Psi.WindowLeftSites + IndexDiff) % PsiLeft.size(); ++i)
               --C;
            
            if (C == PsiLeft.end())
            {
               inplace_qshift(PsiLeft, adjoint(PsiLeft.qshift()));
               C = PsiLeft.begin();
            }
         }
         else if (Index > WindowRightIndex)
         {
            int IndexDiff = Index - WindowRightIndex - 1;

            for (int i = 0; i < (Psi.WindowRightSites + IndexDiff) / PsiRight.size(); ++i)
               inplace_qshift(PsiRight, adjoint(PsiRight.qshift()));

            C = PsiRight.begin();
            for (int i = 0; i < (Psi.WindowRightSites + IndexDiff) % PsiRight.size(); ++i)
               ++C;
         }
         else
         {
            int IndexDiff = Index - WindowLeftIndex;

            C = Psi.Window.begin();
            for (int i = 0; i < IndexDiff; ++i)
               ++C;
         }
      }

      StateComponent operator*() const
      {
         StateComponent Result = *C;

         if (Index == WindowRightIndex + 1)
            Result = Psi.Window.lambda_r() * Psi.Window.RightU() * Result;

         if (Index == WindowLeftIndex)
            Result = Psi.Window.LeftU() * Result;

         return Result;
      }

      ConstIBCIterator& operator++()
      {
         ++Index;
         if (Index < WindowLeftIndex)
         {
            ++C;
            if (C == PsiLeft.end())
            {
               inplace_qshift(PsiLeft, adjoint(PsiLeft.qshift()));
               C = PsiLeft.begin();
            }
         }
         else if (Index == WindowRightIndex + 1)
         {
            C = PsiRight.begin();
            for (int i = 0; i < Psi.WindowRightSites; ++i)
               ++C;
         }
         else if (Index == WindowLeftIndex)
            C = Psi.Window.begin();
         else if (Index <= WindowRightIndex)
            ++C;
         else if (Index > WindowRightIndex + 1)
         {
            ++C;
            if (C == PsiRight.end())
            {
               inplace_qshift(PsiRight, adjoint(PsiRight.qshift()));
               C = PsiRight.begin();
            }
         }

         return *this;
      }

   private:
      IBCWavefunction Psi;
      InfiniteWavefunctionLeft PsiLeft;
      InfiniteWavefunctionRight PsiRight;
      CanonicalWavefunctionBase::const_mps_iterator C;
      int Index;
      int WindowLeftIndex;
      int WindowRightIndex;
};

// NOTE: This function assumes that the left/right boundaries of Psi1 and Psi2
// are identical.
// TODO:
// Move these functions to the correct place.
// Figure out how to properly handle phase.
std::complex<double>
overlap_simple(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2)
{
   int IndexLeft = std::min(Psi1.window_offset(), Psi2.window_offset());
   int IndexRight = std::max(Psi1.window_size() + Psi1.window_offset(),
                             Psi2.window_size() + Psi2.window_offset());

   ConstIBCIterator C1 = ConstIBCIterator(Psi1, IndexLeft);
   ConstIBCIterator C2 = ConstIBCIterator(Psi2, IndexLeft);

   CHECK((*C1).Basis1() == (*C2).Basis1());

   MatrixOperator I = MatrixOperator::make_identity((*C1).Basis1());
   MatrixOperator E = scalar_prod(herm(I), I);

   for (int i = IndexLeft; i <= IndexRight; ++i)
   {
      E = operator_prod(herm(*C1), E, *C2);
      ++C1, ++C2;
   }

   return trace(E);
}

// A more general function to calculate the overlap, which attempts to handle
// the case where the left/right boundaries of Psi1 and Psi2 may be different
// (e.g. if Psi2's boundaries are the complex conjugate of Psi1's).
std::complex<double>
overlap(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2)
{
   CHECK_EQUAL(Psi1.Left.size(), Psi2.Left.size());
   CHECK_EQUAL(Psi1.Right.size(), Psi2.Right.size());

   int IndexLeft1 = Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size());
   int IndexLeft2 = Psi2.window_offset() - ((Psi2.Left.size() - Psi2.WindowLeftSites) % Psi2.Left.size());
   int IndexLeft = std::min(IndexLeft1, IndexLeft2);

   int IndexRight1 = Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size());
   int IndexRight2 = Psi2.window_size() + Psi2.window_offset() + ((Psi2.Right.size() - Psi2.WindowRightSites - 1) % Psi2.Right.size());
   int IndexRight = std::max(IndexRight1, IndexRight2);

   // Before calculating the overlap, we must find the left/right eigenvectors
   // of the left/right boundary transfer matrices.

   // Get the identity quantum number.
   QuantumNumber QI(Psi1.GetSymmetryList());

   // Ensure that the semi-infinite boundaries have the same quantum numbers.
   QuantumNumber QShiftLeft1 = QI;
   for (int i = 0; i < (IndexLeft1 - IndexLeft) / Psi1.Left.size(); ++i)
      QShiftLeft1 = delta_shift(QShiftLeft1, Psi1.Left.qshift());

   QuantumNumber QShiftLeft2 = QI;
   for (int i = 0; i < (IndexLeft2 - IndexLeft) / Psi2.Left.size(); ++i)
      QShiftLeft2 = delta_shift(QShiftLeft2, Psi2.Left.qshift());

   InfiniteWavefunctionLeft Psi1Left = Psi1.Left;
   inplace_qshift(Psi1Left, QShiftLeft1);

   InfiniteWavefunctionLeft Psi2Left = Psi2.Left;
   inplace_qshift(Psi2Left, QShiftLeft2);

   // Calculate the left eigenvector of the left semi-infinite boundary.
   std::complex<double> OverlapL;
   StateComponent ESC;
   std::tie(OverlapL, ESC) = overlap(Psi1Left, Psi2Left, QI);

   // Check that the eigenvalue has magnitude 1.
   //TRACE(OverlapL);
   if (std::abs(std::abs(OverlapL) - 1.0) > OverlapTol)
      WARNING("The overlap of the left boundaries is below threshold.")(OverlapL);

   MatrixOperator E = delta_shift(ESC.front(), adjoint(Psi1.Left.qshift()));

   QuantumNumber QShiftRight1 = QI;
   for (int i = 0; i < (IndexRight - IndexRight1) / Psi1.Right.size(); ++i)
      QShiftRight1 = delta_shift(QShiftRight1, adjoint(Psi1.Right.qshift()));

   QuantumNumber QShiftRight2 = QI;
   for (int i = 0; i < (IndexRight - IndexRight2) / Psi2.Right.size(); ++i)
      QShiftRight2 = delta_shift(QShiftRight2, adjoint(Psi2.Right.qshift()));

   InfiniteWavefunctionRight Psi1Right = Psi1.Right;
   inplace_qshift(Psi1Right, QShiftRight1);

   InfiniteWavefunctionRight Psi2Right = Psi2.Right;
   inplace_qshift(Psi2Right, QShiftRight2);

   // Calculate the right eigenvector of the right semi-infinite boundary.
   std::complex<double> OverlapR;
   StateComponent FSC;
   std::tie(OverlapR, FSC) = overlap(Psi1Right, Psi2Right, QI);

   // Check that the eigenvalue has magnitude 1.
   //TRACE(OverlapR);
   if (std::abs(std::abs(OverlapR) - 1.0) > OverlapTol)
      WARNING("The overlap of the right boundaries is below threshold.")(OverlapR);

   MatrixOperator F = FSC.front();

   // Rescale E and F by the bond dimension.
   E *= std::sqrt(std::min(E.Basis1().total_degree(), E.Basis2().total_degree()));
   F *= std::sqrt(std::min(F.Basis1().total_degree(), F.Basis2().total_degree()));

   if (E.Basis1() == E.Basis2() && F.Basis1() == F.Basis2())
   {
      //E *= std::sqrt(E.Basis1().total_degree());
      //F *= std::sqrt(F.Basis1().total_degree());

      // Remove spurious phase from E and F by setting the phase of the trace
      // to be zero (but only if the trace is nonzero).

      std::complex<double> ETrace = trace(E);

      if (std::abs(ETrace) > TraceTol)
         E *= std::conj(ETrace) / std::abs(ETrace);
      else
         WARNING("The trace of E is below threshold, so the overlap will have a spurious phase contribution.")(ETrace);

      std::complex<double> FTrace = trace(F);

      if (std::abs(FTrace) > TraceTol)
         F *= std::conj(FTrace) / std::abs(FTrace);
      else
         WARNING("The trace of F is below threshold, so the overlap will have a spurious phase contribution.")(FTrace);
   }
   else
      WARNING("Psi1 and Psi2 have different boundary bases, so the overlap will have a spurious phase contribution.");

   // Calculate the overlap.
   ConstIBCIterator C1 = ConstIBCIterator(Psi1, IndexLeft);
   ConstIBCIterator C2 = ConstIBCIterator(Psi2, IndexLeft);

   for (int i = IndexLeft; i <= IndexRight; ++i)
   {
      E = operator_prod(herm(*C1), E, *C2);
      ++C1, ++C2;
   }

   return inner_prod(F, E);
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
      std::string LhsStr, RhsStr;
      int Rotate = 0;
      bool Quiet = false;
      bool Reflect = false;
      bool Conj = false;
      bool Simple = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("cart,c", prog_opt::bool_switch(&ShowCartesian),
          "show the result in cartesian coordinates [equivalent to --real --imag]")
         ("polar,p", prog_opt::bool_switch(&ShowPolar),
          "show the result in polar coodinates [equivalent to --mag --arg]")
         ("real,r", prog_opt::bool_switch(&ShowRealPart),
          "display the real part of the result")
         ("imag,i", prog_opt::bool_switch(&ShowImagPart),
          "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
         ("radians", prog_opt::bool_switch(&ShowRadians),
          "display the argument in radians instead of degrees")
         ("rotate", prog_opt::value(&Rotate),
          "rotate psi2 this many sites to the left before calculating the overlap [default 0]")
         ("reflect", prog_opt::bool_switch(&Reflect),
          "reflect psi2")
         ("conj", prog_opt::bool_switch(&Conj),
          "complex conjugate psi2")
         ("simple", prog_opt::bool_switch(&Simple),
          "assume that psi1 and psi2 have the same semi-infinite boundaries")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the column headings")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&LhsStr), "psi1")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("rhs") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi1> <psi2>\n";
         std::cerr << desc << '\n';
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

      if (Verbose)
         std::cout << "Loading LHS wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi1Ptr;
      if (UseTempFile)
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1Ptr = pheap::ImportHeap(LhsStr);
      }
      else
      {
         Psi1Ptr = pheap::OpenPersistent(LhsStr, mp_pheap::CacheSize(), true);
      }
      IBCWavefunction Psi1 = Psi1Ptr->get<IBCWavefunction>();

      if (Verbose)
         std::cout << "Loading RHS wavefunction...\n";

      pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(RhsStr);
      IBCWavefunction Psi2 = Psi2Ptr->get<IBCWavefunction>();

      if (Verbose)
         std::cout << "Calculating overlap...\n";

      if (Verbose)
         std::cout << "Rotating Psi2 left by " << Rotate << " sites..." << std::endl;
      Psi2.WindowOffset -= Rotate;

      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi2..." << std::endl;
         inplace_reflect(Psi2);
      }

      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi2..." << std::endl;
         inplace_conj(Psi2);
      }

      if (!Quiet)
      {
         int IndexLeft, IndexRight;
         if (Simple)
         {
            IndexLeft = std::min(Psi1.window_offset(), Psi2.window_offset());
            IndexRight = std::max(Psi1.window_size() + Psi1.window_offset(),
                                  Psi2.window_size() + Psi2.window_offset());
         }
         else
         {
            IndexLeft = std::min(Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size()),
                                 Psi2.window_offset() - ((Psi2.Left.size() - Psi2.WindowLeftSites) % Psi2.Left.size()));
            IndexRight = std::max(Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size()),
                                  Psi2.window_size() + Psi2.window_offset() + ((Psi2.Right.size() - Psi2.WindowRightSites - 1) % Psi2.Right.size()));
         }

         std::cout << "#overlap calculated over sites " << IndexLeft << " to " << IndexRight
                   << " (" << IndexRight - IndexLeft + 1 << " sites total)\n";
         if (ShowRealPart)
            std::cout << "#real                   ";
         if (ShowImagPart)
            std::cout << "#imag                   ";
         if (ShowMagnitude)
            std::cout << "#magnitude              ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
         std::cout << '\n';
         std::cout << std::left;
      }

      std::complex<double> x;

      if (Simple)
         x = overlap_simple(Psi1, Psi2);
      else
         x = overlap(Psi1, Psi2);

      PrintFormat(x, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);

      pheap::Shutdown();

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
