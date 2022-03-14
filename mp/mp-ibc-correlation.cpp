// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-correlation.cpp
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
#include "parser/number-parser.h"

// The tolerance for the left/right boundary overlaps for the overlap of two general IBCs.
double const OverlapTol = 1e-12;

#if 0
// The tolerance of the trace of the left/right boundary eigenvectors for
// fixing the phase when calculating the overlap of two general IBCs.
double const TraceTol = 1e-8;
#endif

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

std::complex<double>
overlap(IBCWavefunction const& Psi1, IBCWavefunction const& Psi2,
        MatrixOperator const& E_, MatrixOperator const& F_)
{
   CHECK_EQUAL(Psi1.Left.size(), Psi2.Left.size());
   CHECK_EQUAL(Psi1.Right.size(), Psi2.Right.size());

   int IndexLeft1 = Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size());
   int IndexLeft2 = Psi2.window_offset() - ((Psi2.Left.size() - Psi2.WindowLeftSites) % Psi2.Left.size());
   int IndexLeft = std::min(IndexLeft1, IndexLeft2);

   int IndexRight1 = Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size());
   int IndexRight2 = Psi2.window_size() + Psi2.window_offset() + ((Psi2.Right.size() - Psi2.WindowRightSites - 1) % Psi2.Right.size());
   int IndexRight = std::max(IndexRight1, IndexRight2);

   MatrixOperator E = E_;
   MatrixOperator F = F_;

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
      bool Quiet = false;
      bool Conj = false;
      //bool Simple = false;

      std::string InputPrefix;
      std::string TimestepStr;
      std::string BetastepStr;
      std::string InitialTimeStr;
      std::string InitialBetaStr;
      int N = 0;
      int XMin = 0;
      int XMax = 0;
      int Digits = 0;
      std::complex<double> InitialTime = 0.0;
      // Scaling factor to remove the spurious phase when using the complex conjugate.
      std::complex<double> ScalingFactor = 1.0;

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
         ("conj", prog_opt::bool_switch(&Conj),
          "complex conjugate psi2")
         //("simple", prog_opt::bool_switch(&Simple),
         // "assume that psi1 and psi2 have the same semi-infinite boundaries")
	 ("wavefunction,w", prog_opt::value(&InputPrefix), "Prefix for input wavefunctions")
	 ("timestep,t", prog_opt::value(&TimestepStr), "Timestep")
	 ("num-timesteps,n", prog_opt::value(&N), FormatDefault("Number of timesteps", N).c_str())
	 ("xmin", prog_opt::value(&XMin), FormatDefault("Minimum value of x", XMin).c_str())
	 ("xmax,x", prog_opt::value(&XMax), FormatDefault("Maximum value of x", XMax).c_str())
	 ("precision", prog_opt::value(&Digits), "Decimal precision in time value of the filenames")
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the column headings")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]" << std::endl;
         std::cerr << desc << std::endl;
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

      std::complex<double> Timestep = 0.0;

      if (!TimestepStr.empty())
         Timestep += ParseNumber(TimestepStr);

      Digits = std::max(Digits, formatting::digits(Timestep));

      if (Verbose)
         std::cout << "Loading RHS wavefunction...\n";

      std::string InitialFilename = InputPrefix;

      if (std::real(Timestep) != 0.0)
         InitialFilename += ".t" + formatting::format_digits(InitialTime, Digits);

      if (std::imag(Timestep) != 0.0)
         InitialFilename += ".b" + formatting::format_digits(InitialTime, Digits);

      pvalue_ptr<MPWavefunction> Psi1Ptr;
      if (UseTempFile)
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1Ptr = pheap::ImportHeap(InitialFilename);
      }
      else
      {
         Psi1Ptr = pheap::OpenPersistent(InitialFilename, mp_pheap::CacheSize(), true);
      }
      IBCWavefunction Psi1 = Psi1Ptr->get<IBCWavefunction>();
      IBCWavefunction Psi2;

      MatrixOperator E, F;

      if (Conj)
      {
         Psi2 = Psi1;
         inplace_conj(Psi1);

         int IndexLeft = Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size());

         int IndexRight = Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size());

         // Get the identity quantum number.
         QuantumNumber QI(Psi1.GetSymmetryList());

         // Calculate the left eigenvector of the left semi-infinite boundary.
         std::complex<double> OverlapL;
         StateComponent ESC;
         std::tie(OverlapL, ESC) = overlap(Psi2.Left, Psi1.Left, QI);

         // Check that the eigenvalue has magnitude 1.
         //TRACE(OverlapL);
         if (std::abs(std::abs(OverlapL) - 1.0) > OverlapTol)
            WARNING("The overlap of the left boundaries is below threshold.")(OverlapL);

         E = delta_shift(ESC.front(), adjoint(Psi1.Left.qshift()));

         // Calculate the right eigenvector of the right semi-infinite boundary.
         std::complex<double> OverlapR;
         StateComponent FSC;
         std::tie(OverlapR, FSC) = overlap(Psi2.Right, Psi1.Right, QI);

         // Check that the eigenvalue has magnitude 1.
         //TRACE(OverlapR);
         if (std::abs(std::abs(OverlapR) - 1.0) > OverlapTol)
            WARNING("The overlap of the right boundaries is below threshold.")(OverlapR);

         F = FSC.front();

         // We could try to fix E and F by hand, but it is easier to just force
         // the correlation function for t = 0, x = 0 to be 1.0 by using a
         // scaling factor.
#if 0
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
#endif

         ScalingFactor = 1.0 / overlap(Psi2, Psi1, E, F);
      }

      if (!Quiet)
      {
         if (std::real(Timestep) != 0.0)
            std::cout << "#time         ";
         if (std::imag(Timestep) != 0.0)
            std::cout << "#beta         ";
         if (XMin != 0 || XMax != 0)
            std::cout << "#x            ";
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

      for (int tstep = 0; tstep <= N; ++tstep)
      {
         if (Verbose)
            std::cout << "Loading LHS wavefunction...\n";

         std::string Filename = InputPrefix;
         std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), Digits);
         std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), Digits);

         if (std::real(Timestep) != 0.0)
            Filename += ".t" + TimeStr;

         if (std::imag(Timestep) != 0.0)
            Filename += ".b" + BetaStr;

         pvalue_ptr<MPWavefunction> Psi2Ptr = pheap::ImportHeap(Filename);
         Psi2 = Psi2Ptr->get<IBCWavefunction>();

         Psi2.WindowOffset += (XMin-1);

         for (int x = XMin; x <= XMax; ++x)
         {
            if (Verbose)
               std::cout << "Rotating Psi2 right by " << x << " sites..." << std::endl;
            ++Psi2.WindowOffset;

            if (Verbose)
            {
               int IndexLeft, IndexRight;
               if (!Conj)
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

               std::cout << "Overlap calculated over sites " << IndexLeft << " to " << IndexRight
                         << " (" << IndexRight - IndexLeft + 1 << " sites total)\n";
            }

            std::complex<double> Overlap;

            if (!Conj)
               Overlap = overlap_simple(Psi2, Psi1);
            else
               Overlap = ScalingFactor * overlap(Psi2, Psi1, E, F);

            if (std::real(Timestep) != 0.0)
               std::cout << std::setw(10) << TimeStr << "    ";
            if (std::imag(Timestep) != 0.0)
               std::cout << std::setw(10) << BetaStr << "    ";
            if (XMin != 0 || XMax != 0)
               std::cout << std::setw(10) << x << "    ";

            PrintFormat(Overlap, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
         }
      }
      
      if (Conj)
      {
         std::complex<double> FinalTime = double(N) * Timestep;
         Psi2.WindowOffset -= XMax;

         for (int tstep = 1; tstep <= N; ++tstep)
         {
            if (Verbose)
               std::cout << "Loading RHS wavefunction...\n";

            std::string Filename = InputPrefix;
            std::string TimeStr = formatting::format_digits(std::real(InitialTime + double(tstep)*Timestep), Digits);
            std::string BetaStr = formatting::format_digits(-std::imag(InitialTime + double(tstep)*Timestep), Digits);

            if (std::real(Timestep) != 0.0)
               Filename += ".t" + TimeStr;

            if (std::imag(Timestep) != 0.0)
               Filename += ".b" + BetaStr;

            pvalue_ptr<MPWavefunction> Psi1Ptr = pheap::ImportHeap(Filename);
            IBCWavefunction Psi1 = Psi1Ptr->get<IBCWavefunction>();
            inplace_conj(Psi1);

            TimeStr = formatting::format_digits(std::real(InitialTime + FinalTime + double(tstep)*Timestep), Digits);
            BetaStr = formatting::format_digits(-std::imag(InitialTime + FinalTime + double(tstep)*Timestep), Digits);

            Psi2.WindowOffset += (XMin-1);

            for (int x = XMin; x <= XMax; ++x)
            {

               if (Verbose)
                  std::cout << "Rotating Psi2 right by " << x << " sites..." << std::endl;
               ++Psi2.WindowOffset;

               if (Verbose)
               {
                  int IndexLeft, IndexRight;
                  IndexLeft = std::min(Psi1.window_offset() - ((Psi1.Left.size() - Psi1.WindowLeftSites) % Psi1.Left.size()),
                                       Psi2.window_offset() - ((Psi2.Left.size() - Psi2.WindowLeftSites) % Psi2.Left.size()));
                  IndexRight = std::max(Psi1.window_size() + Psi1.window_offset() + ((Psi1.Right.size() - Psi1.WindowRightSites - 1) % Psi1.Right.size()),
                                        Psi2.window_size() + Psi2.window_offset() + ((Psi2.Right.size() - Psi2.WindowRightSites - 1) % Psi2.Right.size()));

                  std::cout << "Overlap calculated over sites " << IndexLeft << " to " << IndexRight
                            << " (" << IndexRight - IndexLeft + 1 << " sites total)\n";
               }

               std::complex<double> Overlap;

               Overlap = ScalingFactor * overlap(Psi2, Psi1, E, F);

               if (std::real(Timestep) != 0.0)
                  std::cout << std::setw(10) << TimeStr << "    ";
               if (std::imag(Timestep) != 0.0)
                  std::cout << std::setw(10) << BetaStr << "    ";
               if (XMin != 0 || XMax != 0)
                  std::cout << std::setw(10) << x << "    ";

               PrintFormat(Overlap, ShowRealPart, ShowImagPart, ShowMagnitude, ShowArgument, ShowRadians);
            }

            Psi2.WindowOffset -= XMax;
         }
      }

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
