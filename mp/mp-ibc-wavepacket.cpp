// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-ibc-wavepacket.cpp
//
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "common/environment.h"
#include "common/formatting.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/rangelist.h"
#include "common/terminal.h"
#include "common/unique.h"
#include "interface/inittemp.h"
#include "linearalgebra/eigen.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include <boost/version.hpp>
#if BOOST_VERSION >= 107500
#include <boost/math/special_functions/jacobi_theta.hpp>
#endif

namespace prog_opt = boost::program_options;

LinearAlgebra::Vector<std::complex<double>>
operator*(LinearAlgebra::Matrix<std::complex<double>> const& M, LinearAlgebra::Vector<std::complex<double>> const& v)
{
   LinearAlgebra::Vector<std::complex<double>> Result(size1(M));
   for (unsigned i = 0; i < size1(M); ++i)
   {
      std::complex<double> x = 0.0;
      for (unsigned j = 0; j < size2(M); ++j)
      {
         x += M(i,j) * v[j];
      }
      Result[i] = x;
   }
   return Result;
}

// The "wrapped Gaussian" is the sum of the Guassian with mean mu and standard
// deviation sigma and all displacements of this Gaussian by 2*pi.
// This can be calculated in terms of the Jacobi theta function.
double
WrappedGaussian(double x, double mu, double sigma)
{
#if BOOST_VERSION >= 107500
   return 1.0/(2.0*math_const::pi) * boost::math::jacobi_theta3tau(0.5*(x-mu), (0.5*std::pow(sigma,2))/math_const::pi);
#else
   // jacobi_theta.hpp is not in versions of Boost before 1.75, so we perform
   // the sum until machine epsilon is reached, which will not be very
   // expensive for small sigma.

   int JMax = std::ceil(std::sqrt(-2.0*std::log(std::numeric_limits<double>::epsilon()))*sigma);
   if (JMax > 10000)
      WARNING("Attempting to calculate WrappedGaussian with JMax > 10000!");

   double Result = 0;
   for (int j = -JMax; j <= JMax; ++j)
      Result += 1.0/std::sqrt(2.0*math_const::pi)/sigma
         * std::exp(-1.0/2.0/std::pow(sigma,2) * std::pow(std::fmod(x-mu+math_const::pi,2.0*math_const::pi)+(-1.0+2.0*j)*math_const::pi, 2));

   return Result;
#endif
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
   int UCSize = WPVec.size(); // The GS wavefunction unit cell size.
   int WindowSize = WPVec.front().size(); // In terms of unit cells.

   std::vector<LinearWavefunction> PsiWindowVec(WindowSize);

   LinearWavefunction PsiLinearLeft, PsiLinearRight;
   std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);
   std::tie(std::ignore, PsiLinearRight) = get_right_canonical(PsiRight);

   QuantumNumber QShift(PsiLeft.GetSymmetryList());

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

      QShift = delta_shift(QShift, adjoint(PsiLeft.qshift()));

      // Handle all of the middle sites: these are of the form
      // (A B)
      // (0 C)
      for (int i = 1; i < UCSize*WindowSize-1; ++i)
      {
         if (WP == WPCell->end())
         {
            QShift = QuantumNumber(PsiLeft.GetSymmetryList());
            ++WPCell, ++CL, ++CR;
            WP = WPCell->begin();
            PsiWindow = PsiWindowVec.begin();
         }
         StateComponent CLShift = delta_shift(*CL, QShift);
         StateComponent WPShift = delta_shift(*WP, QShift);
         StateComponent CRShift = delta_shift(*CR, QShift);

         StateComponent Z = StateComponent(CLShift.LocalBasis(), CRShift.Basis1(), CLShift.Basis2());
         SumBasis<VectorBasis> NewBasis1(CLShift.Basis2(), WPShift.Basis2());
         SumBasis<VectorBasis> NewBasis2(CLShift.Basis1(), CRShift.Basis1());
         PsiWindow->push_back(tensor_col_sum(tensor_row_sum(CLShift, WPShift, NewBasis1), tensor_row_sum(Z, CRShift, NewBasis1), NewBasis2));

         ++WP, ++PsiWindow;
         QShift = delta_shift(QShift, adjoint(PsiLeft.qshift()));
      }

      // Handle the last site separately: this is of the form
      // (B)
      // (C)
      StateComponent CRShift = delta_shift(*CR, QShift);
      StateComponent WPShift = delta_shift(*WP, QShift);
      SumBasis<VectorBasis> NewBasisBack(WPShift.Basis1(), CRShift.Basis1());
      PsiWindow->push_back(tensor_col_sum(WPShift, CRShift, NewBasisBack));
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
            // calculate the contribution for all x, otherwise, only use the
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
      std::string KStr;
      double KSMA = 0.0;
      std::string KYStr = "0";
      int LatticeUCSize = 1;
      std::string InputPrefix;
      std::string OutputFilename;
      bool Force = false;
      double Sigma = 0.0;
      double KCenter = 0.0;
      double SigmaY = 0.0;
      double KYCenter = 0.0;
      int InputDigits = -1;
      double Tol = 1e-5;
      std::string LeftBoundaryFilename;
      std::string RightBoundaryFilename;
      StatesInfo SInfo;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("wavefunction,w", prog_opt::value(&InputPrefix), "Prefix for input filenames (of the form [prefix].k[k]) [required]")
         ("momentum,k", prog_opt::value(&KStr), "Momentum range of the form start:end:step or start:end,num (in units of pi) [required]")
         ("sma", prog_opt::value(&KSMA), "Use a single-mode approximation using the EA wavefunction for this momentum (in units of pi) [alternative to -k]")
         ("ky", prog_opt::value(&KYStr), "Localize the wavepacker on a cylinder using this y-momentum range (in units of pi)")
         ("digits", prog_opt::value(&InputDigits), "Manually use this number of decimal places for the filenames")
         ("latticeucsize", prog_opt::value(&LatticeUCSize), "Lattice unit cell size [default wavefunction attribute \"LatticeUnitCellSize\" or 1]")
         ("output,o", prog_opt::value(&OutputFilename), "Output filename [required]")
         ("force,f", prog_opt::bool_switch(&Force), "Force overwriting output file")
         ("sigma,s", prog_opt::value(&Sigma), "Convolute with a Gaussian in momentum space with this width (in units of pi)")
         ("kcenter,c", prog_opt::value(&KCenter), FormatDefault("Central momentum of the momentum space Gaussian (in units of pi)", KCenter).c_str())
         ("sigmay", prog_opt::value(&SigmaY), "Convolute with a Gaussian in y-momentum space with this width (in units of pi)")
         ("kycenter", prog_opt::value(&KYCenter), FormatDefault("Central momentum of the y-momentum space Gaussian (in units of pi)", KCenter).c_str())
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance for the wavepacket weight outside the window", Tol).c_str())
         ("min-states", prog_opt::value(&SInfo.MinStates),
          FormatDefault("Minimum number of states to keep", SInfo.MinStates).c_str())
         ("max-states", prog_opt::value<int>(&SInfo.MaxStates),
          FormatDefault("Maximum number of states to keep", SInfo.MaxStates).c_str())
         ("trunc,r", prog_opt::value<double>(&SInfo.TruncationCutoff),
          FormatDefault("Truncation error cutoff", SInfo.TruncationCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&SInfo.EigenvalueCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", SInfo.EigenvalueCutoff).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("wavefunction") == 0 || vm.count("output") == 0 ||
          (vm.count("momentum") == 0) == (vm.count("sma") == 0)) // i.e. momentum XNOR sma
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -w <input-prefix> -o <psi-out>" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      pheap::Initialize(OutputFilename, 1, mp_pheap::PageSize(), mp_pheap::CacheSize(), false, Force);

      RangeList KYList(KYStr);
      RangeList KList;

      if (vm.count("sma") == 0)
      {
         KList = RangeList(KStr);

         if (KList.get_num() <= 1)
         {
            std::cerr << "fatal: The momentum range must contain more than one value." << std::endl;
            return 1;
         }
      }
      else
      {
         KList = RangeList(std::to_string(KSMA));

         // TODO: We might want to be able to use a different SMA for each y-momentum.
         CHECK(vm.count("ky") == 0);

         if (Sigma == 0.0)
         {
            std::cerr << "fatal: --sigma must be specified if using --sma." << std::endl;
            return 1;
         }
      }

      if (InputDigits == -1)
      {
         InputDigits = std::max(formatting::digits(KList.get_start()), formatting::digits(KList.get_step()));
         if (vm.count("kynum"))
            InputDigits = std::max(InputDigits, formatting::digits(KYList.get_step()));
      }

      if (Verbose > 1)
         std::cout << "Loading wavefunctions..." << std::endl;

      // Load first wavefunction: we do this to get the left and right
      // boundaries and the unit cell size.
      InfiniteWavefunctionLeft PsiLeft, PsiLeftOriginal;
      InfiniteWavefunctionRight PsiRight, PsiRightOriginal;
      QuantumNumber LeftQShift, RightQShift;
      int LeftIndex, RightIndex, UCSize;
      {
         std::string InputFilename = InputPrefix;
         if (vm.count("kynum") == 0)
            InputFilename += ".k" + formatting::format_digits(KList.get_start(), InputDigits);
         else
            InputFilename += ".kx" + formatting::format_digits(KList.get_start(), InputDigits)
                           + ".ky" + formatting::format_digits(KYList.get_start(), InputDigits);

         pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
         EAWavefunction Psi = InPsi->get<EAWavefunction>();

         // Get the lattice unit cell size if unspecified.
         if (!vm.count("latticeucsize"))
            LatticeUCSize = InPsi->Attributes()["LatticeUnitCellSize"].get_or_default<int>(1);

         // If the input streams the boundaries, then we save them so we can
         // stream them in the output as well.
         LeftBoundaryFilename = Psi.get_left_filename();
         RightBoundaryFilename = Psi.get_right_filename();

         PsiLeftOriginal = Psi.left();
         PsiRightOriginal = Psi.right();
         LeftQShift = Psi.left_qshift();
         RightQShift = Psi.right_qshift();
         LeftIndex = Psi.left_index();
         RightIndex = Psi.right_index();

         PsiLeft = Psi.left();
         inplace_qshift(PsiLeft, Psi.left_qshift());
         PsiLeft.rotate_left(Psi.left_index());

         PsiRight = Psi.right();
         inplace_qshift(PsiRight, Psi.right_qshift());
         PsiRight.rotate_left(Psi.right_index());

         UCSize = PsiLeft.size();

         CHECK(PsiRight.size() == UCSize);
         CHECK(Psi.window_vec().size() == UCSize);
         CHECK(PsiLeft.qshift() == PsiRight.qshift());
      }

      if (UCSize % LatticeUCSize)
      {
         std::cerr << "fatal: the specified lattice unit cell size must divide the wavefunction unit cell size." << std::endl;
         return 1;
      }

      // The number of lattice unit cells in Psi.
      int LatticeUCsPerPsiUC = UCSize / LatticeUCSize;

      // A vector for each position at the unit cell, which contains a vector
      // of each B matrix for that unit cell position.
      std::vector<std::vector<StateComponent>> BVec(UCSize);
      std::vector<std::complex<double>> ExpIKVec;
      for (double const k : KList)
      {
         for (double const ky : KYList)
         {
            // Load wavefunction.
            std::string InputFilename = InputPrefix;
            if (vm.count("kynum") == 0)
               InputFilename += ".k" + formatting::format_digits(k, InputDigits);
            else
               InputFilename += ".kx" + formatting::format_digits(k, InputDigits)
                              + ".ky" + formatting::format_digits(ky, InputDigits);

            pvalue_ptr<MPWavefunction> InPsi = pheap::ImportHeap(InputFilename);
            EAWavefunction Psi = InPsi->get<EAWavefunction>();

            // We only handle single-site EAWavefunctions at the moment.
            CHECK(Psi.window_size() == 1);

            // TODO: Check each wavefunction has the same left/right boundaries.

            ExpIKVec.push_back(Psi.exp_ik());

            auto BCell = BVec.begin();
            // Here we use the left-gauge fixing condition.
#if 1
            for (WavefunctionSectionLeft Window : Psi.window_vec())
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
            for (WavefunctionSectionLeft Window : Psi.window_vec())
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

      std::vector<std::vector<StateComponent>> WPVec;
      int L; // The final width of the window.

      // If we are using a different B-matrix for each momentum, we need to
      // find the optimal wavepacket numerically.
      if (!vm.count("sma"))
      {
         if (Verbose > 1)
            std::cout << "Localizing wavepacket..." << std::endl;

         // The number of Fourier modes in our momentum space (note that this does
         // not match KNum since we may have missing parts of the spectrum).
         int N = std::round(2.0/KList.get_step()/LatticeUCsPerPsiUC);
         if (std::abs(N*KList.get_step()/2.0 * LatticeUCsPerPsiUC - 1.0) > 0.0)
            std::cerr << "WARNING: Number of Fourier modes " << 2.0/KList.get_step()/LatticeUCsPerPsiUC << " is noninteger! Trying N=" << N << std::endl;

         int NY;
         if (vm.count("kynum"))
         {
            NY = std::round(2.0/KYList.get_step());
            if (std::abs(NY*KYList.get_step()/2.0 - 1.0) > 0.0)
               std::cerr << "WARNING: Number of y Fourier modes " << 2.0/KYList.get_step() << " is noninteger! Trying NY=" << NY << std::endl;
            if (NY < 4)
            {
               std::cerr << "fatal: NY=" << NY << " is less than 4: cannot localize wavepacket along the y-axis!" << std::endl;
               return 1;
            }
         }
         else
            // This makes it such that each Lambda is calculated once.
            NY = 4;

         std::vector<std::complex<double>> FVec;
         std::vector<std::vector<std::vector<std::complex<double>>>> BBVec = CalculateBBVec(BVec);
         // Number of different Fourier modes that we optimize over.
         int Size = KList.get_num()*KYList.get_num();

         int Lambda = 1;
         int LambdaY = 1;
         bool Finished = false;
         while (Lambda < N/2 && !Finished)
         {
            LambdaY = 1;
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
            std::cerr << "fatal: Cannot localize wavepacket." << std::endl;
            return 1;
         }

         if (Verbose > 2)
         {
            // Print the F vector before convolution.
            std::cout << "Printing F before convolution..." << std::endl;
            std::cout << "#kx ky F" << std::endl;
            auto K = KList.begin();
            auto KY = KYList.begin();
            auto F = FVec.begin();
            while (K != KList.end())
            {
               std::cout << *K << " ";
               if (vm.count("ky"))
                  std::cout << *KY << " ";
               std::cout << formatting::format_complex(*F) << std::endl;
               ++K, ++KY, ++F;
            }
         }

         // Convolute with momentum space Gaussian.
         if (Sigma != 0.0)
         {
            if (Verbose > 1)
               std::cout << "Convoluting with momentum space Gaussian..." << std::endl;
            auto K = KList.begin();
            auto F = FVec.begin();
            while (K != KList.end())
            {
               // Since we have a periodic domain in momentum space, we cannot
               // just use a normal Gaussian, so we use the "wrapped Gaussian",
               // which is defined on a periodic domain.
               *F *= WrappedGaussian(math_const::pi*(*K), math_const::pi*KCenter, math_const::pi*Sigma);
               ++K, ++F;
            }
         }

         // Convolute with y-momentum space Gaussian.
         if (SigmaY != 0.0)
         {
            if (Verbose > 1)
               std::cout << "Convoluting with y-momentum space Gaussian..." << std::endl;
            auto KY = KYList.begin();
            auto F = FVec.begin();
            while (KY != KYList.end())
            {
               *F *= WrappedGaussian(math_const::pi*(*KY), math_const::pi*KYCenter, math_const::pi*SigmaY);
               ++KY, ++F;
            }
         }

         if (Verbose > 2)
         {
            // Print the F vector after convolution.
            std::cout << "Printing F after convolution..." << std::endl;
            std::cout << "#kx ky F" << std::endl;
            auto K = KList.begin();
            auto KY = KYList.begin();
            auto F = FVec.begin();
            while (K != KList.end())
            {
               std::cout << *K << " ";
               if (vm.count("ky"))
                  std::cout << *KY << " ";
               std::cout << formatting::format_complex(*F) << std::endl;
               ++K, ++KY, ++F;
            }
         }

         if (Verbose > 1)
         {
            std::vector<std::vector<StateComponent>> WPVecFull = CalculateWPVec(BVec, ExpIKVec, FVec, N/2);
            // Print the norm of the B-matrices in WPVec for each unit cell.
            std::cout << "Printing wavepacket B-matrix norms..." << std::endl;
            std::cout << "#i j norm_frob(B(j)[i])" << std::endl;
            int i = 0;
            for (auto const& WPCell : WPVecFull)
            {
               int j = -N/2;
               for (auto const& WP : WPCell)
                  std::cout << i << " " << j++ << " " << norm_frob(WP) << std::endl;
               ++i;
            }
         }

         // Calculate the new value of Lambda post-convolution.
         int LambdaNew;
         LinearAlgebra::Vector<std::complex<double>> FVector(FVec.begin(), FVec.end());
         for (LambdaNew = Lambda; LambdaNew < N/2; ++LambdaNew)
         {
            LinearAlgebra::Matrix<std::complex<double>> NLambdaMat
               = CalculateNLambda(BBVec, ExpIKVec, N, LambdaNew, LambdaY, LatticeUCSize);
            double Error = std::real(inner_prod(FVector, NLambdaMat * FVector));

            if (Verbose > 2)
               std::cout << "LambdaNew = " << LambdaNew
                         << ", Error = " << Error << std::endl;

            if (Error < Tol)
               break;
         }

         if (Verbose > 1)
            std::cout << "Updating Lambda to " << LambdaNew << std::endl;

         Lambda = LambdaNew;

         if (Verbose > 1)
            std::cout << "Constructing window..." << std::endl;

         WPVec = CalculateWPVec(BVec, ExpIKVec, FVec, Lambda);
         L = Lambda;
      }
      // If we are using the single-mode approximation, the B-matrix is the
      // same for each momentum, and we can localize the wavepacket exactly.
      else
      {
         // Scale Sigma and KCenter by pi and the number of lattice unit cells in Psi.
         double SigmaScale = math_const::pi * LatticeUCsPerPsiUC * Sigma;
         double KCenterScale = math_const::pi * LatticeUCsPerPsiUC * KCenter;

         // The number of Psi unit cells in the window.
         L = std::ceil(std::sqrt(-2.0*std::log(Tol))/SigmaScale);

         if (Verbose > 1)
            std::cout << "Constructing window with L = " << L << std::endl;

         WPVec = std::vector<std::vector<StateComponent>>(UCSize);

         for (int n = -L; n <= L; ++n)
         {
            auto B = BVec.begin();
            auto WP = WPVec.begin();
            while (B != BVec.end())
            {
               WP->push_back(std::exp(-0.5*std::pow(SigmaScale*n,2)) * std::exp(std::complex<double>(0.0,KCenterScale*n)) * B->front());
               ++B, ++WP;
            }
         }
      }

      LinearWavefunction PsiWindowLinear = ConstructPsiWindow(PsiLeft, PsiRight, WPVec);

      if (Verbose > 1)
         std::cout << "Orthogonalizing window..." << std::endl;

      MatrixOperator I = MatrixOperator::make_identity(PsiWindowLinear.Basis2());
      WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), I, Verbose-1);

      if (Verbose > 1)
         std::cout << "Normalizing and truncating window..." << std::endl;

      MatrixOperator LambdaWindow;
      std::tie(PsiWindowLinear, LambdaWindow) = get_left_canonical(PsiWindow);
      LambdaWindow *= 1.0 / norm_frob(LambdaWindow);

      PsiWindowLinear.set_back(prod(PsiWindowLinear.get_back(), LambdaWindow));

      truncate_left_orthogonal(PsiWindowLinear, SInfo, Verbose-1);

      PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(PsiWindowLinear), I, Verbose-1);

      if (Verbose > 1)
         std::cout << "Saving wavefunction..." << std::endl;

      LeftQShift = delta_shift(LeftQShift, PsiLeft.qshift());

      for (int i = 0; i < WPVec.front().size(); ++i)
         RightQShift = delta_shift(RightQShift, adjoint(PsiRight.qshift()));

      IBCWavefunction PsiOut(PsiLeftOriginal, PsiWindow, PsiRightOriginal, LeftQShift, RightQShift, -L*UCSize,
                             (PsiLeft.size() - LeftIndex) % PsiLeft.size(), RightIndex);

      // Stream the boundaries, if the input files do.
      PsiOut.set_left_filename(LeftBoundaryFilename);
      PsiOut.set_right_filename(RightBoundaryFilename);

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
