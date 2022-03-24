// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-excitation-ansatz.cpp
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

#define LEFT_RIGHT 1

#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "lattice/latticesite.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/lanczos.h"
#include "lattice/infinitelattice.h"
#include "lattice/infinite-parser.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "wavefunction/operator_actions.h"
#include "mp-algorithms/gmres.h"

#define EXP_I(k) exp(std::complex<double>(0.0, 1.0) * (k))

namespace prog_opt = boost::program_options;

double
norm_frob_sq(std::deque<MatrixOperator> const& Input)
{
   double Result = 0.0;
   for (auto I : Input)
      Result += norm_frob_sq(I);
   return Result;
}

double
norm_frob(std::deque<MatrixOperator> const& Input)
{
   return std::sqrt(norm_frob_sq(Input));
}

std::deque<MatrixOperator>&
operator*=(std::deque<MatrixOperator>& Input, double x)
{
   for (auto& I : Input)
      I *= x;
   return Input;
}

std::deque<MatrixOperator>
operator*(double x, std::deque<MatrixOperator> const& Input)
{
   std::deque<MatrixOperator> Result = Input;
   Result *= x;
   return Result;
}

std::deque<MatrixOperator>&
operator+=(std::deque<MatrixOperator>& Input1, std::deque<MatrixOperator> const& Input2)
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

std::deque<MatrixOperator>&
operator-=(std::deque<MatrixOperator>& Input1, std::deque<MatrixOperator> const& Input2)
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

std::complex<double>
inner_prod(std::deque<MatrixOperator> const& Input1, std::deque<MatrixOperator> const& Input2)
{
   std::complex<double> Result = 0.0;
   auto I1 = Input1.begin();
   auto I2 = Input2.begin();
   while (I1 != Input1.end())
   {
      Result += inner_prod(*I1, *I2);
      ++I1, ++I2;
   }
   return Result;
}

// NB: This is a modified version of mp-algorithms/lanczos.h which returns the
// found eigenvalues, and not just the lowest one.
// TODO: Use a better version of Lanczos to find the lowest n eigenvalues.
template <typename VectorType, typename MultiplyFunctor>
LinearAlgebra::Vector<double>
LanczosFull(VectorType& Guess, MultiplyFunctor MatVecMultiply, int& Iterations,
        double& Tol, int MinIter = 2, int Verbose = 0)
{
   std::vector<VectorType>     v;          // the Krylov vectors
   std::vector<VectorType>     Hv;         // the Krylov vectors

   LinearAlgebra::Matrix<double> SubH(Iterations+1, Iterations+1, 0.0);

   VectorType w = Guess;

   double Beta = norm_frob(w);
   //TRACE(Beta);
   CHECK(!std::isnan(Beta));
   // double OrigBeta = Beta;      // original norm, not used
   w *= 1.0 / Beta;
   v.push_back(w);

   w = MatVecMultiply(v[0]);
   Hv.push_back(w);
   SubH(0,0) = real(inner_prod(v[0], w));
   w -= SubH(0,0) * v[0];

   Beta = norm_frob(w);
   //TRACE(Beta);
   if (Beta < LanczosBetaTol)
   {
      if (Verbose > 0)
         std::cerr << "lanczos: immediate return, invariant subspace found, Beta="
                   << Beta << '\n';
      Guess = v[0];
      Iterations = 1;
      Tol = Beta;
      LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,1),
                                             LinearAlgebra::range(0,1));
      LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
      return EValues;
   }

   // It isn't meaningful to do a Krylov algorithm with only one matrix-vector multiply,
   // but we postpone the check until here to allow the corner case of
   // Iterations==1 but the algorithm converged in 1 step,
   // which might happen for example if the Hilbert space is 1-dimensional
   // (external checks that set Iterations=max(Iterations,Dimension) are not unreasonable).
   CHECK(Iterations > 1)("Number of iterations must be greater than one")(Iterations);

   for (int i = 1; i < Iterations; ++i)
   {
      SubH(i, i-1) = SubH(i-1, i) = Beta;
      w *= 1.0 / Beta;
      v.push_back(w);
      w = MatVecMultiply(v[i]);
      Hv.push_back(w);
      w -= Beta*v[i-1];
      SubH(i,i) = real(inner_prod(v[i], w));
      w -= SubH(i,i) * v[i];
      Beta = norm_frob(w);
      //TRACE(Beta);


      if (Beta < LanczosBetaTol)
      {
         // Early return, we can't improve over the previous energy and eigenvector
         if (Verbose > 0)
            std::cerr << "lanczos: early return, invariant subspace found, Beta="
                      << Beta << ", iterations=" << (i+1) << '\n';
         Iterations = i+1;
         LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i+1),
                                                LinearAlgebra::range(0,i+1));
         LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
         double Theta = EValues[0];    // smallest eigenvalue
         double SpectralDiameter = EValues[i] - EValues[0];  // largest - smallest
         VectorType y = M(0,0) * v[0];
         for (int j = 1; j <= i; ++j)
            y += M(0,j) * v[j];
         Tol = Beta / SpectralDiameter;
         Guess = y;
         return EValues;
      }

      // solution of the tridiagonal subproblem
      LinearAlgebra::Matrix<double> M = SubH(LinearAlgebra::range(0,i+1),
                                             LinearAlgebra::range(0,i+1));
#if 0
      if (std::isnan(M(0,0)))
      {
         std::ofstream Out("lanczos_debug.txt");
         Out << "NAN encountered in Lanczos\n"
             << "Beta=" << Beta << "\n\n"
             << "norm_frob(Guess)=" << norm_frob(Guess) << "\n\n"
             << "Guess=" << Guess << "\n\n"
             << "M=" << "\n\n"
             << "SubH=" << SubH << "\n\n";
         for (unsigned n = 0; n < v.size(); ++n)
         {
            Out << "V[" << n << "]=" << v[n] << "\n\n";
         }
      }
#endif

      LinearAlgebra::Vector<double> EValues = DiagonalizeHermitian(M);
      //for (double E : EValues)
      //   std::cout << E << std::endl;
      double Theta = EValues[0];    // smallest eigenvalue
      std::cout << Theta << std::endl;
      double SpectralDiameter = EValues[i] - EValues[0];  // largest - smallest
      VectorType y = M(0,0) * v[0];
      for (int j = 1; j <= i; ++j)
         y += M(0,j) * v[j];

      // residual = H*y - Theta*y
      VectorType r = (-Theta) * y;
      for (int j = 0; j < i; ++j)
         r += M(0,j) * Hv[j];

      double ResidNorm = norm_frob(r);

      if (ResidNorm < fabs(Tol * SpectralDiameter) && i+1 >= MinIter)
         //if (ResidNorm < Tol && i+1 >= MinIter)
      {
         if (Verbose > 0)
            std::cerr << "lanczos: early return, residual norm below threshold, ResidNorm="
                      << ResidNorm << ", SpectralDiameter=" << SpectralDiameter
                      << ", iterations=" << (i+1) << '\n';
         Iterations = i+1;
         Tol = ResidNorm / SpectralDiameter;
         Guess = y;
         return EValues;
      }
      else if (Verbose > 2)
      {
         std::cerr << "lanczos: Eigen=" << Theta << ", ResidNorm="
                      << ResidNorm << ", SpectralDiameter=" << SpectralDiameter
                      << ", iterations=" << (i+1) << '\n';
      }

      if (i == Iterations-1) // finished?
      {
         Guess = y;
         Tol = -ResidNorm / SpectralDiameter;
         return EValues;
      }
   }

   PANIC("Should never get here");
   LinearAlgebra::Vector<double> EValues(0);
   return EValues;
}


struct HEff
{
#if LEFT_RIGHT
   HEff(InfiniteWavefunctionLeft const& PsiLeft_, InfiniteWavefunctionLeft const& PsiRight_,
        BasicTriangularMPO const& HamMPO_, double k_, double GMRESTol_, int Verbose_)
      : PsiLeft(PsiLeft_), PsiRight(PsiRight_),
        HamMPO(HamMPO_), k(k_), GMRESTol(GMRESTol_), Verbose(Verbose_)
   {
      CHECK_EQUAL(PsiLeft.size(), PsiRight.size());

      std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);

      std::tie(U, D, PsiLinearRight) = get_right_canonical(PsiRight);

      PsiLinearRight.set_front(prod(U, PsiLinearRight.get_front()));

      for (StateComponent C : PsiLinearLeft)
         NullLeftDeque.push_back(NullSpace2(C));

      if (HamMPO.size() < PsiLeft.size())
         HamMPO = repeat(HamMPO, PsiLeft.size() / HamMPO.size());

      BlockHamL = Initial_E(HamMPO, PsiLeft.Basis1());
      std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPO, GMRESTol, Verbose);
      if (Verbose > 0)
         std::cout << "Left energy = " << LeftEnergy << std::endl;

      //BlockHamL.back() -= LeftEnergy * BlockHamL.front();

      BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

      BlockHamR = Initial_F(HamMPO, PsiLinearRight.Basis2());
      MatrixOperator Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
      Rho = delta_shift(Rho, adjoint(PsiRight.qshift()));

      std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiLinearRight, PsiRight.qshift(), HamMPO, Rho, GMRESTol, Verbose);
      if (Verbose > 0)
         std::cout << "Right energy = " << RightEnergy << std::endl;

      //BlockHamR.front() -= RightEnergy * BlockHamR.back();

      BlockHamR = delta_shift(BlockHamR, PsiRight.qshift());

      TRACE(inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamR));

      BlockHamR.front() -= (RightEnergy + inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamR)) * BlockHamR.back();

      TRACE(inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamR));

      BlockHamLDeque.push_back(BlockHamL);
      auto CL = PsiLinearLeft.begin();
      auto O = HamMPO.begin();
      while (CL != PsiLinearLeft.end())
      {
         BlockHamLDeque.push_back(contract_from_left(*O, herm(*CL), BlockHamLDeque.back(), *CL));
         ++CL, ++O;
      }

      BlockHamRDeque.push_front(BlockHamR);
      auto CR = PsiLinearRight.end();
      O = HamMPO.end();
      while (CR != PsiLinearRight.begin())
      {
         --CR, --O;
         BlockHamRDeque.push_front(contract_from_right(herm(*O), *CR, BlockHamRDeque.front(), herm(*CR)));
      }

      TRACE(inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamLDeque.back(), PsiLeft.lambda_r())), BlockHamRDeque.back()));
   }

   std::deque<MatrixOperator>
   operator()(std::deque<MatrixOperator> const& XDeque) const
   {
      std::deque<StateComponent> BDeque;
      auto NL = NullLeftDeque.begin();
      auto X = XDeque.begin();
      while (NL != NullLeftDeque.end())
      {
         BDeque.push_back(prod(*NL, *X));
         ++NL, ++X;
      }
      
      StateComponent BL, BR;
      
      //MatrixOperator Rho = U * D * adjoint(U);
      MatrixOperator Rho = D;

#if 0
      TRACE(norm_frob(Rho - inject_left(Rho, PsiLinearRight, PsiLinearLeft)));
      TRACE(norm_frob(Rho - inject_right(Rho, PsiLinearRight, PsiLinearLeft)));
      TRACE(norm_frob(Rho - inject_left(Rho, PsiLinearLeft, PsiLinearRight)));
      TRACE(norm_frob(Rho - inject_right(Rho, PsiLinearLeft, PsiLinearRight)));
#endif

      std::complex<double> LeftEnergy = SolveSimpleMPO_Left2(BL, BlockHamL, PsiLinearLeft, PsiLinearRight*EXP_I(k), BDeque,
                                                             PsiLeft.qshift(), HamMPO, Rho, Rho, GMRESTol, Verbose);
      std::complex<double> RightEnergy = SolveSimpleMPO_Right2(BR, BlockHamR, PsiLinearLeft, PsiLinearRight*EXP_I(k), BDeque,
                                                               PsiRight.qshift(), HamMPO, Rho, Rho, GMRESTol, Verbose);

      //TRACE(norm_frob(BL))(norm_frob(BR));

      std::deque<MatrixOperator> Result;
      auto B = BDeque.begin();
      NL = NullLeftDeque.begin();
      auto O = HamMPO.begin();
      auto BHL = BlockHamLDeque.begin();
      auto BHR = BlockHamRDeque.begin();
      ++BHR;
      while (B != BDeque.end())
      {
         Result.push_back(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *NL)), *BHR));
         //TRACE(norm_frob(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *NL)), *BHR)));
         TRACE(trace(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *B)), *BHR)));
         ++B, ++NL, ++O, ++BHL, ++BHR;
      }

      StateComponent Tmp = BL;
      B = BDeque.begin();
      NL = NullLeftDeque.begin();
      auto CL = PsiLinearLeft.begin();
      auto CR = PsiLinearRight.begin();
      O = HamMPO.begin();
      BHL = BlockHamLDeque.begin();
      BHR = BlockHamRDeque.begin();
      ++BHR;
      auto R = Result.begin();
      while (B != BDeque.end())
      {
         *R += scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *NL)), *BHR);
         //TRACE(norm_frob(scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *NL)), *BHR)));
         TRACE(trace(scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *B)), *BHR)));
         Tmp = contract_from_left(*O, herm(*CR), Tmp, *CL) + contract_from_left(*O, herm(*B), *BHL, *CL);
         ++B, ++NL, ++CL, ++CR, ++O, ++BHL, ++BHR, ++R;
      }

      Tmp = BR;
      B = BDeque.end();
      NL = NullLeftDeque.end();
      CL = PsiLinearLeft.end();
      CR = PsiLinearRight.end();
      O = HamMPO.end();
      BHL = BlockHamLDeque.end();
      --BHL;
      BHR = BlockHamRDeque.end();
      R = Result.end();
      X = XDeque.end();
      while (B != BDeque.begin())
      {
         --B, --NL, --CL, --CR, --O, --BHL, --BHR, --R;
         --X;
         *R += scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *NL)), Tmp);
         //TRACE(norm_frob(scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *NL)), Tmp)));
         TRACE(trace(scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *B)), Tmp)));
         TRACE(inner_prod(*R, *X));
         Tmp = contract_from_right(herm(*O), *CL, Tmp, herm(*CR)) + contract_from_right(herm(*O), *B, *BHR, herm(*CR));
      }

      //TRACE(norm_frob(Result));
      return Result;
   }
#else
   HEff(InfiniteWavefunctionLeft const& PsiLeft_, InfiniteWavefunctionLeft const& PsiRight_,
        BasicTriangularMPO const& HamMPO_, double k_, double GMRESTol_, int Verbose_)
      : PsiLeft(PsiLeft_), PsiRight(PsiRight_),
        HamMPO(HamMPO_), k(k_), GMRESTol(GMRESTol_), Verbose(Verbose_)
   {
      CHECK_EQUAL(PsiLeft.size(), PsiRight.size());

      std::tie(PsiLinearLeft, D) = get_left_canonical(PsiLeft);
      std::tie(PsiLinearRight, std::ignore) = get_left_canonical(PsiRight);

      for (StateComponent C : PsiLinearLeft)
         NullLeftDeque.push_back(NullSpace2(C));

      if (HamMPO.size() < PsiLeft.size())
         HamMPO = repeat(HamMPO, PsiLeft.size() / HamMPO.size());

      BlockHamL = Initial_E(HamMPO, PsiLeft.Basis1());
      std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPO, GMRESTol, Verbose);
      if (Verbose > 0)
         std::cout << "Left energy = " << LeftEnergy << std::endl;

      //BlockHamL.back() -= LeftEnergy * BlockHamL.front();

      BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

      BlockHamR = Initial_F(HamMPO, PsiLinearLeft.Basis2());
      MatrixOperator Ident = BlockHamR.back();

      MatrixOperator Rho = D*D;
      Rho = delta_shift(Rho, adjoint(PsiLeft.qshift()));
      BlockHamR.back() = Rho;

      std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiLinearLeft, PsiLeft.qshift(), HamMPO, Ident, GMRESTol, Verbose);
      if (Verbose > 0)
         std::cout << "Right energy = " << RightEnergy << std::endl;

      //BlockHamR.front() -= RightEnergy * BlockHamR.back();

      BlockHamR = delta_shift(BlockHamR, PsiRight.qshift());

      TRACE(inner_prod(BlockHamL, BlockHamR));

      BlockHamR.front() -= (RightEnergy + inner_prod(BlockHamL, BlockHamR)) * BlockHamR.back();

      TRACE(inner_prod(BlockHamL, BlockHamR));

      BlockHamLDeque.push_back(BlockHamL);
      auto CL = PsiLinearLeft.begin();
      auto O = HamMPO.begin();
      while (CL != PsiLinearLeft.end())
      {
         BlockHamLDeque.push_back(contract_from_left(*O, herm(*CL), BlockHamLDeque.back(), *CL));
         ++CL, ++O;
      }

      BlockHamRDeque.push_front(BlockHamR);
      auto CR = PsiLinearLeft.end();
      O = HamMPO.end();
      while (CR != PsiLinearLeft.begin())
      {
         --CR, --O;
         BlockHamRDeque.push_front(contract_from_right(herm(*O), *CR, BlockHamRDeque.front(), herm(*CR)));
      }

      TRACE(inner_prod(BlockHamLDeque.back(), BlockHamRDeque.back()));
   }

   std::deque<MatrixOperator>
   operator()(std::deque<MatrixOperator> const& XDeque) const
   {
      std::deque<StateComponent> BDeque;
      auto NL = NullLeftDeque.begin();
      auto X = XDeque.begin();
      while (NL != NullLeftDeque.end())
      {
         BDeque.push_back(prod(*NL, *X));
         ++NL, ++X;
      }
      
      StateComponent BL, BR;
      
      MatrixOperator Ident = Initial_F(HamMPO, PsiLinearLeft.Basis2()).back();
      MatrixOperator Rho = D*D;

      std::complex<double> LeftEnergy = SolveSimpleMPO_Left2(BL, BlockHamL, PsiLinearLeft, PsiLinearLeft*EXP_I(k), BDeque,
                                                             PsiLeft.qshift(), HamMPO, Ident, Rho, GMRESTol, Verbose);
      std::complex<double> RightEnergy = SolveSimpleMPO_Right2(BR, BlockHamR, PsiLinearLeft, PsiLinearLeft*EXP_I(k), BDeque,
                                                               PsiLeft.qshift(), HamMPO, Ident, Rho, GMRESTol, Verbose);

      //TRACE(norm_frob(BL))(norm_frob(BR));

      std::deque<MatrixOperator> Result;
      auto B = BDeque.begin();
      NL = NullLeftDeque.begin();
      auto O = HamMPO.begin();
      auto BHL = BlockHamLDeque.begin();
      auto BHR = BlockHamRDeque.begin();
      ++BHR;
      while (B != BDeque.end())
      {
         Result.push_back(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *NL)), *BHR));
         //TRACE(norm_frob(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *NL)), *BHR)));
         TRACE(trace(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *B)), *BHR)));
         ++B, ++NL, ++O, ++BHL, ++BHR;
      }

      StateComponent Tmp = BL;
      B = BDeque.begin();
      NL = NullLeftDeque.begin();
      auto CL = PsiLinearLeft.begin();
      auto CR = PsiLinearLeft.begin();
      O = HamMPO.begin();
      BHL = BlockHamLDeque.begin();
      BHR = BlockHamRDeque.begin();
      ++BHR;
      auto R = Result.begin();
      while (B != BDeque.end())
      {
         *R += scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *NL)), *BHR);
         //TRACE(norm_frob(scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *NL)), *BHR)));
         TRACE(trace(scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *B)), *BHR)));
         Tmp = contract_from_left(*O, herm(*CR), Tmp, *CL) + contract_from_left(*O, herm(*B), *BHL, *CL);
         ++B, ++NL, ++CL, ++CR, ++O, ++BHL, ++BHR, ++R;
      }

      Tmp = BR;
      B = BDeque.end();
      NL = NullLeftDeque.end();
      CL = PsiLinearLeft.end();
      CR = PsiLinearLeft.end();
      O = HamMPO.end();
      BHL = BlockHamLDeque.end();
      --BHL;
      BHR = BlockHamRDeque.end();
      R = Result.end();
      X = XDeque.end();
      while (B != BDeque.begin())
      {
         --B, --NL, --CL, --CR, --O, --BHL, --BHR, --R;
         --X;
         *R += scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *NL)), Tmp);
         //TRACE(norm_frob(scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *NL)), Tmp)));
         TRACE(trace(scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *B)), Tmp)));
         TRACE(inner_prod(*R, *X));
         Tmp = contract_from_right(herm(*O), *CL, Tmp, herm(*CR)) + contract_from_right(herm(*O), *B, *BHR, herm(*CR));
      }

      //TRACE(norm_frob(Result));
      return Result;
   }
#endif

   std::deque<MatrixOperator>
   InitialGuess()
   {
      std::deque<MatrixOperator> Result;
      auto NL = NullLeftDeque.begin();
      auto CR = PsiLinearRight.begin();
      while (NL != NullLeftDeque.end())
      {
         MatrixOperator C = MakeRandomMatrixOperator((*NL).Basis2(), (*CR).Basis2());
         C *= 1.0 / norm_frob(C);
         Result.push_back(C);
         ++NL, ++CR;
      }

      return Result;
   }

   InfiniteWavefunctionLeft const& PsiLeft;
   InfiniteWavefunctionLeft const& PsiRight;
   BasicTriangularMPO HamMPO;
   double k;
   double GMRESTol;
   int Verbose;
   LinearWavefunction PsiLinearLeft, PsiLinearRight;
   StateComponent BlockHamL, BlockHamR;
   std::deque<StateComponent> BlockHamLDeque, BlockHamRDeque;
   std::deque<StateComponent> NullLeftDeque;
   MatrixOperator U;
   RealDiagonalOperator D;
};

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string InputFileLeft;
      std::string InputFileRight;
      std::string HamStr;
      double k = 0;
      bool Random = false;
      bool Force = false;
      double GMRESTol = 1E-13;    // tolerance for GMRES for the initial H matrix elements.
      int Iter = 50;
      int MinIter = 4;
      double Tol = 1E-16;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("left,l", prog_opt::value(&InputFileLeft),
          "Input iMPS wavefunction for the left semi-infinite strip [required]")
         ("right,r", prog_opt::value(&InputFileRight),
          "Input iMPS wavefunction for the right semi-infinite strip")
         ("Hamiltonian,H", prog_opt::value(&HamStr),
          "Operator to use for the Hamiltonian")
         ("momentum,k", prog_opt::value(&k),
          "Excitation momentum")
         ("maxiter", prog_opt::value<int>(&Iter),
          FormatDefault("Maximum number of Lanczos iterations", Iter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations", MinIter).c_str())
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Error tolerance for the Lanczos eigensolver", Tol).c_str())
         ("gmrestol", prog_opt::value(&GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", GMRESTol).c_str())
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("left") == 0 || vm.count("Hamiltonian") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      mp_pheap::InitializeTempPHeap();

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(InputFileLeft);
      pvalue_ptr<MPWavefunction> InPsiRight;

      if (vm.count("right"))
         InPsiRight = pheap::ImportHeap(InputFileRight);
      else
         InPsiRight = pheap::ImportHeap(InputFileLeft);

      InfiniteWavefunctionLeft PsiLeft = InPsiLeft->get<InfiniteWavefunctionLeft>();

      // There are some situations where the first method does not work
      // properly, so temporarily use the second method as a workaround.
#if 0
      InfiniteWavefunctionRight PsiRight = InPsiRight->get<InfiniteWavefunctionLeft>();
#else
      InfiniteWavefunctionLeft PsiRightLeft = InPsiRight->get<InfiniteWavefunctionLeft>();

      MatrixOperator U;
      RealDiagonalOperator D;
      LinearWavefunction PsiRightLinear;
      std::tie(U, D, PsiRightLinear) = get_right_canonical(PsiRightLeft);

      InfiniteWavefunctionRight PsiRight(U*D, PsiRightLinear, PsiRightLeft.qshift());
#endif

      // If the bases of the two boundary unit cells have only one quantum
      // number sector, manually ensure that they match.
      // FIXME: This workaround probably will not work for non-Abelian symmetries.
      if (PsiLeft.Basis2().size() == 1 && PsiRight.Basis1().size() == 1)
         inplace_qshift(PsiRight, delta_shift(PsiLeft.Basis2()[0], adjoint(PsiRight.Basis1()[0])));

      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO, HamMPOLeft, HamMPORight;
      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

      HEff EffectiveHamiltonian(PsiLeft, PsiRightLeft, HamMPO, k, GMRESTol, Verbose);

      std::deque<MatrixOperator> XDeque = EffectiveHamiltonian.InitialGuess();

      LinearAlgebra::Vector<double> EValues = LanczosFull(XDeque, EffectiveHamiltonian, Iter, Tol, MinIter, Verbose);

      std::cout << "Eigenvalues:" << std::endl;

      for (double E : EValues)
         std::cout << E << std::endl;
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
