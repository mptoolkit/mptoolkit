// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-ibc-wavepacket.cpp
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
#include "lattice/infinite-parser.h"
#include "lattice/infinitelattice.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"

namespace prog_opt = boost::program_options;

std::vector<std::vector<StateComponent>>
CalculateWPVec(std::vector<std::vector<StateComponent>> const& BVec, std::vector<std::complex<double>> const& ExpIKVec,
               std::vector<std::complex<double>> const& FVec, int const N)
{
   // A vector for each position at the unit cell, which contains a vector
   // of each wavepacket B matrix for that unit cell position.
   std::vector<std::vector<StateComponent>> WPVec(BVec.size());

   auto WPCell = WPVec.begin();
   auto BCell = BVec.begin();
   while (BCell != BVec.end())
   {
      auto F = FVec.begin();
      auto ExpIK = ExpIKVec.begin();
      auto B = BCell->begin();
      *WPCell = std::vector<StateComponent>(N, StateComponent(B->LocalBasis(), B->Basis1(), B->Basis2()));
      while (B != BCell->end())
      {
         // Leftmost unit cell position of the N-unit cell window.
         int j = -N/2;
         for (auto& WP : *WPCell)
         {
            WP += *F * std::pow(*ExpIK, j) * *B;
            ++j;
         }
         ++F, ++ExpIK, ++B;
      }
      ++WPCell, ++BCell;
   }

   return WPVec;
}

std::vector<std::complex<double>>
CalculateNLambdaF(std::vector<std::vector<StateComponent>> const& BVec, std::vector<std::complex<double>> const& ExpIKVec,
                  std::vector<std::complex<double>> const& FVec, int const N, int const Lambda)
{
   std::vector<std::complex<double>> NLambdaFVec(FVec.size(), 0.0);

   auto BCell = BVec.begin();
   while (BCell != BVec.end())
   {
      auto FJ = FVec.begin();
      auto ExpIKJ = ExpIKVec.begin();
      auto BJ = BCell->begin();
      while (BJ != BCell->end())
      {
         auto NLambdaFI = NLambdaFVec.begin();
         auto ExpIKI = ExpIKVec.begin();
         auto BI = BCell->begin();
         while (BI != BCell->end())
         {
            std::complex<double> Coeff = 0.0;
            for (int j = Lambda+1; j < N-Lambda; ++j)
               Coeff += std::pow(std::conj(*ExpIKI) * *ExpIKJ, j);

            *NLambdaFI += *FJ * inner_prod(*BI, *BJ) * Coeff;
            ++NLambdaFI, ++ExpIKI, ++BI;
         }
         ++FJ, ++ExpIKJ, ++BJ;
      }
      ++BCell;
   }

   return NLambdaFVec;
}

struct NLambdaFunctor
{
   NLambdaFunctor(std::vector<std::vector<StateComponent>> const& BVec_,
                  std::vector<std::complex<double>> const& ExpIKVec_,
                  int const N_, int const Lambda_)
      : BVec(BVec_), ExpIKVec(ExpIKVec_), N(N_), Lambda(Lambda_)
   {
   }

   std::vector<std::complex<double>> operator()(std::vector<std::complex<double>> const& FVec) const
   {
      return CalculateNLambdaF(BVec, ExpIKVec, FVec, N, Lambda);
   }

   void SetLambda(int const Lambda_)
   {
      Lambda = Lambda_;
   }

   std::vector<std::vector<StateComponent>> const& BVec;
   std::vector<std::complex<double>> const& ExpIKVec;
   int N;
   int Lambda;
};

std::vector<std::complex<double>>&
operator*=(std::vector<std::complex<double>>& Input, std::complex<double> x)
{
   for (auto& I : Input)
      I *= x;
   return Input;
}

std::vector<std::complex<double>>
operator*(std::complex<double> x, std::vector<std::complex<double>> const& Input)
{
   std::vector<std::complex<double>> Result = Input;
   Result *= x;
   return Result;
}

std::vector<std::complex<double>>&
operator+=(std::vector<std::complex<double>>& Input1, std::vector<std::complex<double>> const& Input2)
{
   auto I1 = Input1.begin();
   auto I2 = Input2.begin();
   while (I1 != Input1.end())
   {
      *I1 += *I2;
      ++I1, ++I2;
   }
   return Input1;
}

std::vector<std::complex<double>>
operator+(std::vector<std::complex<double>> const& Input1, std::vector<std::complex<double>> const& Input2)
{
   std::vector<std::complex<double>> Result = Input1;
   Result += Input2;
   return Result;
}

std::vector<std::complex<double>>&
operator-=(std::vector<std::complex<double>>& Input1, std::vector<std::complex<double>> const& Input2)
{
   auto I1 = Input1.begin();
   auto I2 = Input2.begin();
   while (I1 != Input1.end())
   {
      *I1 -= *I2;
      ++I1, ++I2;
   }
   return Input1;
}

std::vector<std::complex<double>>
operator-(std::vector<std::complex<double>> const& Input1, std::vector<std::complex<double>> const& Input2)
{
   std::vector<std::complex<double>> Result = Input1;
   Result -= Input2;
   return Result;
}

double
norm_frob_sq(std::vector<std::complex<double>> const& Input)
{
   double Result = 0.0;
   for (auto I : Input)
      Result += std::norm(I);
   return Result;
}

double
norm_frob(std::vector<std::complex<double>> const& Input)
{
   return std::sqrt(norm_frob_sq(Input));
}

std::complex<double>
inner_prod(std::vector<std::complex<double>> const& Input1, std::vector<std::complex<double>> const& Input2)
{
   std::complex<double> Result = 0.0;
   auto I1 = Input1.begin();
   auto I2 = Input2.begin();
   while (I1 != Input1.end())
   {
      Result += std::conj(*I1) * *I2;
      ++I1, ++I2;
   }
   return Result;
}

struct IdentityFunctor
{
   IdentityFunctor()
   {
   }

   std::vector<std::complex<double>> operator()(std::vector<std::complex<double>> const& FVec) const
   {
      return FVec;
   }
};

struct InnerProdFunctor
{
   InnerProdFunctor()
   {
   }

   std::complex<double> operator()(std::vector<std::complex<double>> const& FVec1, std::vector<std::complex<double>> const& FVec2) const
   {
      return inner_prod(FVec1, FVec2);
   }
};

// FIXME: We have to include this header here so that it uses the functions
// defined above. There must be a better solution.
#include "mp-algorithms/conjugategradient.h"

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double KMax = 0;
      double KMin = 0;
      int KNum = 1;
      std::string InputPrefix;
      int InputDigits = 0;
      std::string OutputFilename;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("kmax,k", prog_opt::value(&KMax), FormatDefault("Maximum momentum (divided by pi)", KMax).c_str())
         ("kmin", prog_opt::value(&KMin), FormatDefault("Minimum momentum (divided by pi)", KMin).c_str())
         ("knum", prog_opt::value(&KNum), "Number of momentum steps to use")
         ("wavefunction,w", prog_opt::value(&InputPrefix), "Prefix for input filenames (of the form [prefix].k[k])")
         ("output,o", prog_opt::value(&OutputFilename), "Output filename")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help"))
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      mp_pheap::InitializeTempPHeap();

      double KStep = (KMax-KMin)/(KNum-1);
      InputDigits = std::max(formatting::digits(KMax), formatting::digits(KStep));

      // Read input wavefunctions.
      std::vector<EAWavefunction> PsiVec;
      for (int n = 0; n < KNum; ++n)
      {
         std::string InputFilename = InputPrefix + ".k" + formatting::format_digits(KMin + KStep*n, InputDigits);
         pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
         PsiVec.push_back(InPsi->get<EAWavefunction>());
         // We only handle single-site EAWavefunctions at the moment.
         CHECK(PsiVec.back().window_size() == 1);
      }

      // A vector for each position at the unit cell, which contains a vector
      // of each B matrix for that unit cell position.
      std::vector<std::vector<StateComponent>> BSymVec(PsiVec.back().Left.size());
      std::vector<std::complex<double>> ExpIKVec;

      for (EAWavefunction Psi : PsiVec)
      {
         ExpIKVec.push_back(Psi.ExpIK);

         auto BSym = BSymVec.begin();
         // Here we use the left-gauge fixing condition.
#if 1
         for (WavefunctionSectionLeft Window : Psi.WindowVec)
         {
            LinearWavefunction PsiLinear;
            MatrixOperator U;
            std::tie(PsiLinear, U) = get_left_canonical(Window);
            // Note that we assume that the window is single-site.
            BSym->push_back(PsiLinear.get_front()*U);
            ++BSym;
         }
         // TODO: Symmetric gauge-fixing.
#else
         // Get the right null space matrices corresponding to each A-matrix in PsiRight.
         LinearWavefunction PsiLinearLeft;
         std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(Psi.Left);

         LinearWavefunction PsiLinearRight;
         std::tie(std::ignore, PsiLinearRight) = get_right_canonical(Psi.Right);

         std::vector<StateComponent> NullRightVec;
         for (StateComponent C : PsiLinearRight)
            NullRightVec.push_back(NullSpace1(C));

         auto NR = NullRightVec.begin();
         auto AL = PsiLinearLeft.begin();
         auto AR = PsiLinearRight.begin();
         for (WavefunctionSectionLeft Window : Psi.WindowVec)
         {
            LinearWavefunction PsiLinear;
            MatrixOperator U;
            std::tie(PsiLinear, U) = get_left_canonical(Window);
            // Note that we assume that the window is single-site.
            StateComponent BL = PsiLinear.get_front()*U;

            // Find the B-matrix satisfying the right-gauge fixing condition.
            MatrixOperator XR = scalar_prod(BL, herm(*NR));
            StateComponent BR = prod(XR, *NR);
            // Scale norm to match BL
            //BR *= norm_frob(BL) / norm_frob(BR);

            TRACE(inner_prod(BR, *AR))(inner_prod(BR, *AL));
            TRACE(inner_prod(BL, *AL))(inner_prod(BL, *AR));
            TRACE(inner_prod(BR, BL));
            TRACE(norm_frob(BL))(norm_frob(BR))(norm_frob(XR));

            BSym->push_back(0.5*(BL+BR));
            ++BSym;
            ++NR, ++AR, ++AL;
         }
#endif
      }

      // The number of Fourier modes in our momentum space (note that this does
      // not match KNum since we may have missing parts of the spectrum).
      // FIXME: This way of calculating N is bound to cause problems.
      int N = std::round(2.0/KStep);
      TRACE(N)(2.0/KStep);

      int Lambda = 5;

      std::vector<std::complex<double>> FVec(ExpIKVec.size(), 1.0);
      NLambdaFunctor NLambda(BSymVec, ExpIKVec, N, Lambda);
      int MaxIter = 100;
      double Tol = 1e-5;

      TRACE(norm_frob(FVec));
      TRACE(inner_prod(FVec, NLambda(FVec)));
      ConjugateGradient(FVec, NLambda, std::vector<std::complex<double>>(FVec.size(), 0.0), MaxIter, Tol, IdentityFunctor(), InnerProdFunctor(), InnerProdFunctor());
      TRACE(MaxIter)(Tol)(norm_frob(FVec));
      FVec *= std::sqrt(ExpIKVec.size()) / norm_frob(FVec);
      TRACE(inner_prod(FVec, NLambda(FVec)));

      std::vector<std::vector<StateComponent>> WPVec = CalculateWPVec(BSymVec, ExpIKVec, FVec, 2*N);
      for (auto const& WPCell : WPVec)
      {
         int j = -N;
         for (auto const& WP : WPCell)
            std::cout << j++ << " " << norm_frob(WP) << std::endl;
      }
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
