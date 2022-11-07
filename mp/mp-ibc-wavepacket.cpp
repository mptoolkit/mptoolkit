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
#include "linearalgebra/eigen.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include <boost/math/special_functions/jacobi_theta.hpp>

namespace prog_opt = boost::program_options;

// The "wrapped Gaussian" is the sum of the Guassian with mean mu and standard
// deviation sigma and all displacements of this Gaussian by 2*pi.
// This can be calculated in terms of the Jacobi theta function.
double
WrappedGaussian(double x, double mu, double sigma)
{
   return 1.0/(2.0*math_const::pi)*boost::math::jacobi_theta3tau(0.5*(math_const::pi*(x-mu)),(0.5*std::pow(sigma,2))/math_const::pi);
}

// Calculate the B-matrices for a wavepacket using sampling function FVec.
std::vector<std::vector<StateComponent>>
CalculateWPVec(std::vector<std::vector<StateComponent>> const& BVec, std::vector<std::complex<double>> const& ExpIKVec,
               std::vector<std::complex<double>> const& FVec, int const Lambda)
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
      *WPCell = std::vector<StateComponent>(2*Lambda+1, StateComponent(B->LocalBasis(), B->Basis1(), B->Basis2()));
      while (B != BCell->end())
      {
         // Leftmost unit cell position of the N-unit cell window.
         int j = -Lambda;
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

// Construct the triangular MPS window for the wavepacket.
LinearWavefunction
ConstructPsiWindow(InfiniteWavefunctionLeft PsiLeft, InfiniteWavefunctionRight PsiRight,
                   std::vector<std::vector<StateComponent>> const& WPVec)
{
   int UCSize = WPVec.front().size(); // The GS wavefunction unit cell size.
   int WindowSize = WPVec.size(); // In terms of unit cells.

   std::vector<LinearWavefunction> PsiWindowVec(UCSize);

   LinearWavefunction PsiLinearLeft, PsiLinearRight;
   std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);
   std::tie(std::ignore, PsiLinearRight) = get_right_canonical(PsiRight);

   // Just in case we somehow have a single site window.
   if (UCSize == 1 && WindowSize == 1)
      PsiWindowVec.front().push_back(WPVec.front().front());
   else
   {
      auto CL = PsiLinearLeft.begin();
      auto CR = PsiLinearRight.begin();
      auto WPCell = WPVec.begin();
      auto WP = WPCell->begin();
      auto PsiWindow = PsiWindowVec.begin();

      // Handle the first site separately: this is of the form
      // (A B)
      SumBasis<VectorBasis> NewBasisFront((*CL).Basis2(), (*WP).Basis2());
      PsiWindow->push_back(tensor_row_sum(*CL, *WP, NewBasisFront));
      ++WP, ++PsiWindow;

      // Handle all of the middle sites: these are of the form
      // (A B)
      // (0 C)
      for (int i = 1; i < UCSize*WindowSize-1; ++i)
      {
         if (WP == WPCell->end())
         {
            ++WPCell, ++CL, ++CR;
            WP = WPCell->begin();
            PsiWindow = PsiWindowVec.begin();
         }
         StateComponent Z = StateComponent((*CL).LocalBasis(), (*CR).Basis1(), (*CL).Basis2());
         SumBasis<VectorBasis> NewBasis1((*CL).Basis2(), (*WP).Basis2());
         SumBasis<VectorBasis> NewBasis2((*CL).Basis1(), (*CR).Basis1());
         PsiWindow->push_back(tensor_col_sum(tensor_row_sum(*CL, *WP, NewBasis1), tensor_row_sum(Z, *CR, NewBasis1), NewBasis2));
         ++WP, ++PsiWindow;
      }

      // Handle the last site separately: this is of the form
      // (B)
      // (C)
      SumBasis<VectorBasis> NewBasisBack((*WP).Basis1(), (*CR).Basis1());
      PsiWindow->push_back(tensor_col_sum(*WP, *CR, NewBasisBack));
   }

   // Concatenate PsiWindowVec.
   LinearWavefunction PsiWindowVecFull;
   for (auto const& PsiWindow : PsiWindowVec)
      PsiWindowVecFull.push_back(PsiWindow);

   return PsiWindowVecFull;
}

// Calculate the "N_Lambda" matrix from Eq. (A6) in Van Damme et al., Phys. Rev. Research 3, 013078.
LinearAlgebra::Matrix<std::complex<double>>
CalculateNLambda(std::vector<std::vector<std::vector<std::complex<double>>>> const& BBVec, std::vector<std::complex<double>> const& ExpIKVec,
                  int const N, int const Lambda, int const LambdaY, int const LatticeUCSize)
{
   int Size = ExpIKVec.size();
   LinearAlgebra::Matrix<std::complex<double>> NLambdaMat(Size, Size, 0.0);

   auto BBCell = BBVec.begin();
   for (int m = 0; m < BBVec.size(); ++m)
   {
      auto ExpIKJ = ExpIKVec.begin();
      auto BBJ = BBCell->begin();
      for (int Col = 0; Col < Size; ++Col)
      {
         auto ExpIKI = ExpIKVec.begin();
         auto BBI = BBJ->begin();
         for (int Row = 0; Row < Size; ++Row)
         {
            std::complex<double> Coeff = 0.0;
            // If LambdaY < m % LatticeUCSize < LatticeUCSize - LambdaY,
            // calculate the contriution for all x, otherwise, only use the
            // contribution from Lambda < x < N - Lambda.
            if ((m % LatticeUCSize) > LambdaY && (m % LatticeUCSize) < LatticeUCSize - LambdaY)
               for (int j = 0; j < N; ++j)
                  Coeff += std::pow(std::conj(*ExpIKI) * *ExpIKJ, j);
            else
               for (int j = Lambda+1; j < N-Lambda; ++j)
                  Coeff += std::pow(std::conj(*ExpIKI) * *ExpIKJ, j);

            NLambdaMat(Col, Row) += *BBI * Coeff;
            ++ExpIKI, ++BBI;
         }
         ++ExpIKJ, ++BBJ;
      }
      ++BBCell;
   }

   return NLambdaMat;
}

// Calculate the vector containing matrices for the inner products <B_i, B_j>:
// since this matrix only needs to be calculated once, doing it separately
// makes the solver slightly faster.
std::vector<std::vector<std::vector<std::complex<double>>>>
CalculateBBVec(std::vector<std::vector<StateComponent>> const& BVec)
{
   std::vector<std::vector<std::vector<std::complex<double>>>> BBVec(BVec.size(),
      std::vector<std::vector<std::complex<double>>>(BVec.front().size(),
      std::vector<std::complex<double>>(BVec.front().size())));

   auto BCell = BVec.begin();
   auto BBCell = BBVec.begin();
   while (BCell != BVec.end())
   {
      auto BJ = BCell->begin();
      auto BBJ = BBCell->begin();
      while (BJ != BCell->end())
      {
         auto BI = BCell->begin();
         auto BBI = BBJ->begin();
         while (BI != BCell->end())
         {
            *BBI = inner_prod(*BI, *BJ);
            ++BI, ++BBI;
         }
         ++BJ, ++BBJ;
      }
      ++BCell, ++BBCell;
   }

   return BBVec;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      double KMax = 1.0;
      double KMin = 0.0;
      int KNum = 1;
      double KYMax = 1.0;
      double KYMin = 0.0;
      int KYNum = 1;
      int LatticeUCSize = 1;
      std::string InputPrefix;
      std::string OutputFilename;
      double Sigma = 0.0;
      double KCenter = 0.0;
      double SigmaY = 0.0;
      double KYCenter = 0.0;
      int InputDigits = -1;
      double Tol = 1e-5;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("kmax", prog_opt::value(&KMax), FormatDefault("Maximum momentum (divided by pi)", KMax).c_str())
         ("kmin", prog_opt::value(&KMin), FormatDefault("Minimum momentum (divided by pi)", KMin).c_str())
         ("knum", prog_opt::value(&KNum), "Number of momentum steps to use [required]")
         ("kymax", prog_opt::value(&KYMax), FormatDefault("Maximum y-momentum (divided by pi)", KYMax).c_str())
         ("kymin", prog_opt::value(&KYMin), FormatDefault("Minimum y-momentum (divided by pi)", KYMin).c_str())
         ("kynum", prog_opt::value(&KYNum), "Number of y-momentum steps to use")
         ("latticeucsize", prog_opt::value(&LatticeUCSize), "Size of lattice unit cell")
         ("wavefunction,w", prog_opt::value(&InputPrefix), "Prefix for input filenames (of the form [prefix].k[k]) [required]")
         ("output,o", prog_opt::value(&OutputFilename), "Output filename [required]")
         ("digits", prog_opt::value(&InputDigits), "Manually use this number of decimal places for the filenames")
         ("sigma,s", prog_opt::value(&Sigma), "Convolute with a Gaussian in momentum space with this width")
         ("kcenter,k", prog_opt::value(&KCenter), FormatDefault("Central momentum of the momentum space Gaussian", KCenter).c_str())
         ("sigmay", prog_opt::value(&SigmaY), "Convolute with a Gaussian in y-momentum space with this width")
         ("kycenter", prog_opt::value(&KYCenter), FormatDefault("Central momentum of the y-momentum space Gaussian", KCenter).c_str())
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance for the wavepacket weight outside the [-Lambda, Lambda] window", Tol).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("output") == 0 || vm.count("knum") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      mp_pheap::InitializeTempPHeap();

      CHECK(KNum > 1);

      double KStep = (KMax-KMin)/(KNum-1);
      double KYStep = KYNum == 1 ? 0.0 : (KYMax-KYMin)/(KYNum-1);

      if (InputDigits == -1)
      {
         InputDigits = std::max(formatting::digits(KMax), formatting::digits(KStep));
         if (vm.count("kynum"))
            InputDigits = std::max(InputDigits, formatting::digits(KYStep));
      }

      if (Verbose > 1)
         std::cout << "Loading wavefunctions..." << std::endl;

      // Load first wavefunction: we do this to get the left and right
      // boundaries and the unit cell size.
      InfiniteWavefunctionLeft PsiLeft;
      InfiniteWavefunctionRight PsiRight;
      int UCSize;
      {
         std::string InputFilename = InputPrefix;
         if (vm.count("kynum") == 0)
            InputFilename += ".k" + formatting::format_digits(KMin, InputDigits);
         else
            InputFilename += ".kx" + formatting::format_digits(KMin, InputDigits)
                           + ".ky" + formatting::format_digits(KYMin, InputDigits);

         pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
         EAWavefunction Psi = InPsi->get<EAWavefunction>();

         PsiLeft = Psi.Left;
         PsiRight = Psi.Right;
         UCSize = PsiLeft.size();

         CHECK(PsiRight.size() == UCSize);
         CHECK(Psi.WindowVec.size() == UCSize);
      }

      // A vector for each position at the unit cell, which contains a vector
      // of each B matrix for that unit cell position.
      std::vector<std::vector<StateComponent>> BVec(UCSize);
      std::vector<std::complex<double>> ExpIKVec;
      std::vector<double> KVec;
      std::vector<double> KYVec;
      for (int n = 0; n < KNum; ++n)
      {
         for (int m = 0; m < KYNum; ++m)
         {
            // Load wavefunction.
            std::string InputFilename = InputPrefix;
            if (vm.count("kynum") == 0)
               InputFilename += ".k" + formatting::format_digits(KMin + KStep*n, InputDigits);
            else
               InputFilename += ".kx" + formatting::format_digits(KMin + KStep*n, InputDigits)
                              + ".ky" + formatting::format_digits(KYMin + KYStep*m, InputDigits);

            pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
            EAWavefunction Psi = InPsi->get<EAWavefunction>();

            // We only handle single-site EAWavefunctions at the moment.
            CHECK(Psi.window_size() == 1);

            // TODO: Check each wavefunction has the same left/right boundaries.

            ExpIKVec.push_back(Psi.ExpIK);
            KVec.push_back(KMin + KStep*n);
            KYVec.push_back(KYMin + KYStep*m);

            auto BCell = BVec.begin();
            // Here we use the left-gauge fixing condition.
#if 1
            for (WavefunctionSectionLeft Window : Psi.WindowVec)
            {
               LinearWavefunction PsiLinear;
               MatrixOperator U;
               std::tie(PsiLinear, U) = get_left_canonical(Window);
               // Note that we assume that the window is single-site.
               BCell->push_back(PsiLinear.get_front()*U);
               ++BCell;
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

               BCell->push_back(0.5*(BL+BR));
               ++BCell;
               ++NR, ++AR, ++AL;
            }
#endif
         }
      }

      if (Verbose > 1)
         std::cout << "Localizing wavepacket..." << std::endl;

      // The number of Fourier modes in our momentum space (note that this does
      // not match KNum since we may have missing parts of the spectrum).
      int N = std::round(2.0/KStep);
      if (N*KStep/2.0 - 1.0 > std::numeric_limits<double>::epsilon())
         std::cerr << "WARNING: 2/KStep=" << 2.0/KStep << " is noninteger! Trying N=" << N << std::endl;

      int NY;
      if (vm.count("kynum"))
      {
         NY = std::round(2.0/KYStep);
         if (NY*KYStep/2.0 - 1.0 > std::numeric_limits<double>::epsilon())
            std::cerr << "WARNING: 2/KYStep=" << 2.0/KYStep << " is noninteger! Trying NY=" << NY << std::endl;
         if (NY < 4)
         {
            std::cerr << "FATAL: NY=" << NY << " is less than 4: cannot localize wavepacket along the y-axis." << std::endl;
            return 1;
         }
      }
      else
         // This makes it such that each Lambda is calculated once.
         NY = 4;


      std::vector<std::complex<double>> FVec;
      std::vector<std::vector<std::vector<std::complex<double>>>> BBVec = CalculateBBVec(BVec);
      // Number of different Fourier modes that we optimize over.
      int Size = KNum*KYNum;

      int Lambda = 1;
      bool Finished = false;
      while (Lambda < N/2 && !Finished)
      {
         int LambdaY = 1;
         // Only try values of LambdaY up to the current value of Lambda.
         // If we aren't localising along the y-axis, then this loop will only
         // run once, since we set NY/2 = 2.
         while (LambdaY < std::min(Lambda+1, NY/2) && !Finished)
         {
            LinearAlgebra::Matrix<std::complex<double>> NLambdaMat = CalculateNLambda(BBVec, ExpIKVec, N, Lambda, LambdaY, LatticeUCSize);
            LinearAlgebra::Vector<double> EValues = LinearAlgebra::DiagonalizeHermitian(NLambdaMat);

            if (Verbose > 0)
            {
               std::cout << "Lambda=" << Lambda;
               if (vm.count("kynum"))
                  std::cout << " LambdaY=" << LambdaY;
               std::cout << " ExternalWeight=" << std::real(EValues[0])
                         << std::endl;
            }

            if (std::real(EValues[0]) < Tol)
            {
               Finished = true;
               // Extract the smallest eigenvector.
               FVec = std::vector<std::complex<double>>(NLambdaMat.data(), NLambdaMat.data()+Size);
            }

            ++LambdaY;
         }
         ++Lambda;
      }

      if (!Finished)
      {
         std::cerr << "FATAL: Cannot localize wavepacket." << std::endl;
         return 1;
      }

      // Convolute with momentum space Gaussian.
      if (Sigma != 0.0)
      {
         if (Verbose > 1)
            std::cout << "Convoluting with momentum space Gaussian..." << std::endl;
         auto K = KVec.begin();
         auto F = FVec.begin();
         while (K != KVec.end())
         {
            //*F *= std::exp(1/2.0/std::pow(Sigma, 2)*(std::pow(math_const::pi*(*K-KCenter),2)));
            // Since we have a periodic domain in momentum space, we cannot
            // just use a normal Gaussian, so we use the "wrapped Gaussian",
            // which is defined on a periodic domain.
            *F *= WrappedGaussian(*K, KCenter, Sigma);
            ++K, ++F;
         }
      }

      // Convolute with y-momentum space Gaussian.
      if (SigmaY != 0.0)
      {
         if (Verbose > 1)
            std::cout << "Convoluting with y-momentum space Gaussian..." << std::endl;
         auto KY = KYVec.begin();
         auto F = FVec.begin();
         while (KY != KYVec.end())
         {
            *F *= WrappedGaussian(*KY, KYCenter, SigmaY);
            ++KY, ++F;
         }
      }

      if (Verbose > 1)
      {
         // Print the norm of the B-matrices in WPVec for each unit cell.
         std::vector<std::vector<StateComponent>> WPVec = CalculateWPVec(BVec, ExpIKVec, FVec, N/2);
         std::cout << "Printing wavepacket B-matrix norms..." << std::endl;
         std::cout << "#i j norm_frob(B(j)[i])" << std::endl;
         int i = 0;
         for (auto const& WPCell : WPVec)
         {
            int j = -N/2;
            for (auto const& WP : WPCell)
               std::cout << i << " " << j++ << " " << norm_frob(WP) << std::endl;
            ++i;
         }
      }

      if (Verbose > 1)
         std::cout << "Constructing window..." << std::endl;

      std::vector<std::vector<StateComponent>> WPVec = CalculateWPVec(BVec, ExpIKVec, FVec, Lambda);
      LinearWavefunction PsiWindowLinear = ConstructPsiWindow(PsiLeft, PsiRight, WPVec);

      if (Verbose > 1)
         std::cout << "Orthogonalizing window..." << std::endl;

      MatrixOperator I = MatrixOperator::make_identity(PsiWindowLinear.Basis2());
      WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), I, Verbose-1);

      if (Verbose > 1)
         std::cout << "Normalizing window..." << std::endl;

      MatrixOperator LambdaWindow;
      std::tie(PsiWindowLinear, LambdaWindow) = get_left_canonical(PsiWindow);
      LambdaWindow *= 1.0 / norm_frob(LambdaWindow);
      PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), LambdaWindow, Verbose-1);

      if (Verbose > 1)
         std::cout << "Saving wavefunction..." << std::endl;

      IBCWavefunction PsiOut(PsiLeft, PsiWindow, PsiRight, -Lambda*UCSize);

      MPWavefunction Wavefunction;
      Wavefunction.Wavefunction() = std::move(PsiOut);
      Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
      Wavefunction.SetDefaultAttributes();
      pvalue_ptr<MPWavefunction> PsiPtr(new MPWavefunction(Wavefunction));
      pheap::ExportHeap(OutputFilename, PsiPtr);
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
