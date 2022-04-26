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
#include "common/unique.h"

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
      double Theta = EValues[0];    // smallest eigenvalue

      if (Verbose > 0)
         std::cout << "i=" << i
                   << ", E0=" << Theta << std::endl;

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
   HEff() {}

   // Initializer for the case where the A-matrices to the left and the right
   // of the excitation correspond to the same ground state Psi.
   HEff(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& HamMPO_,
        double k, double GMRESTol_, int Verbose_)
      : PsiLeft(Psi_), PsiRight(Psi_),
        HamMPO(HamMPO_), GMRESTol(GMRESTol_), Verbose(Verbose_)
   {
      ExpIK = exp(std::complex<double>(0.0, math_const::pi) * k);

      // Get PsiLeft and PsiRight as LinearWavefunctions.
      std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);

      MatrixOperator U;
      RealDiagonalOperator D;
      std::tie(U, D, PsiLinearRight) = get_right_canonical(PsiRight);
      PsiLinearRight.set_front(prod(U, PsiLinearRight.get_front()));

      // Get the leading eigenvectors for the mixed transfer matrix of PsiLeft
      // and PsiRight: for use with SolveSimpleMPOLeft/Right2.
      RhoL = delta_shift(PsiLeft.lambda_r(), PsiLeft.qshift());
      RhoR = PsiLeft.lambda_r();

      // Solve the left Hamiltonian environment.
      BlockHamL = Initial_E(HamMPO, PsiLeft.Basis1());
      std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPO, GMRESTol, Verbose-1);
      if (Verbose > 0)
         std::cout << "Left energy = " << LeftEnergy << std::endl;

      BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

      // Solve the right Hamiltonian environment.
      BlockHamR = Initial_F(HamMPO, PsiLinearRight.Basis2());
      MatrixOperator Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
      Rho = delta_shift(Rho, adjoint(PsiRight.qshift()));

      std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiLinearRight, PsiRight.qshift(),
                                                              HamMPO, Rho, GMRESTol, Verbose-1);
      if (Verbose > 0)
         std::cout << "Right energy = " << RightEnergy << std::endl;

      // Remove the contribution from the ground state energy density.
      BlockHamR.front() -= (RightEnergy + inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamR)) * BlockHamR.back();

      this->Initialize();
   }

   // Initializer for the case where PsiLeft and PsiRight are two DIFFERENT ground states.
   HEff(InfiniteWavefunctionLeft const& PsiLeft_, InfiniteWavefunctionLeft const& PsiRight_,
        BasicTriangularMPO const& HamMPO_, double k, double GMRESTol_, int Verbose_)
      : PsiLeft(PsiLeft_), PsiRight(PsiRight_),
        HamMPO(HamMPO_), GMRESTol(GMRESTol_), Verbose(Verbose_)
   {
      CHECK_EQUAL(PsiLeft.size(), PsiRight.size());
      CHECK_EQUAL(PsiLeft.qshift(), PsiRight.qshift());

      ExpIK = exp(std::complex<double>(0.0, math_const::pi) * k);

      // Get PsiLeft and PsiRight as LinearWavefunctions.
      std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);

      MatrixOperator U;
      RealDiagonalOperator D;
      std::tie(U, D, PsiLinearRight) = get_right_canonical(PsiRight);
      PsiLinearRight.set_front(prod(U, PsiLinearRight.get_front()));

      // Since the leading eigenvalue of the left/right mixed transfer matrix
      // has magnitude < 1, we do not need to orthogonalize the E/F matrix
      // elements against its eigenvectors.
      RhoL = MatrixOperator();
      RhoR = MatrixOperator();

      // Solve the left Hamiltonian environment.
      BlockHamL = Initial_E(HamMPO, PsiLeft.Basis1());
      std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPO, GMRESTol, Verbose-1);
      if (Verbose > 0)
         std::cout << "Left energy = " << LeftEnergy << std::endl;

      BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

      // Solve the right Hamiltonian environment.
      BlockHamR = Initial_F(HamMPO, PsiLinearRight.Basis2());
      MatrixOperator Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
      Rho = delta_shift(Rho, adjoint(PsiRight.qshift()));

      std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiLinearRight, PsiRight.qshift(),
                                                              HamMPO, Rho, GMRESTol, Verbose-1);
      if (Verbose > 0)
         std::cout << "Right energy = " << RightEnergy << std::endl;

      // TODO: Figure out how to shift the overall constant in BlockHamR such
      // that the energy corresponds to the excitation energy above the ground
      // state.
#if 0
      // Solve the right Hamiltonian environment for PsiLeft (to fix the total energy).
      LinearWavefunction PsiLinear;
      std::tie(U, D, PsiLinear) = get_right_canonical(PsiLeft);
      PsiLinear.set_front(prod(U, PsiLinear.get_front()));

      StateComponent BlockHamLR = Initial_F(HamMPO, PsiLinear.Basis2());
      Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
      Rho = delta_shift(Rho, adjoint(PsiLeft.qshift()));

      SolveSimpleMPO_Right(BlockHamLR, PsiLinear, PsiLeft.qshift(), HamMPO, Rho, GMRESTol, Verbose-1);

      // Remove the contribution from the ground state energy density.
      BlockHamR.front() -= (RightEnergy + inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamLR)) * BlockHamR.back();
#endif
      BlockHamR.front() -= 2.0 * RightEnergy * BlockHamR.back();

      this->Initialize();
   }

   // This function performs the part of the initialization common to both cases above.
   void
   Initialize()
   {
      // Get the null space matrices corresponding to each A-matrix in PsiLeft.
      for (StateComponent C : PsiLinearLeft)
         NullLeftDeque.push_back(NullSpace2(C));

      if (HamMPO.size() < PsiLeft.size())
         HamMPO = repeat(HamMPO, PsiLeft.size() / HamMPO.size());

      // Construct the partially contracted left Hamiltonian environments in the unit cell.
      BlockHamLDeque.push_back(BlockHamL);
      auto CL = PsiLinearLeft.begin();
      auto O = HamMPO.begin();
      while (CL != PsiLinearLeft.end())
      {
         BlockHamLDeque.push_back(contract_from_left(*O, herm(*CL), BlockHamLDeque.back(), *CL));
         ++CL, ++O;
      }

      // Same for the right environments.
      BlockHamRDeque.push_front(BlockHamR);
      auto CR = PsiLinearRight.end();
      O = HamMPO.end();
      while (CR != PsiLinearRight.begin())
      {
         --CR, --O;
         BlockHamRDeque.push_front(contract_from_right(herm(*O), *CR, BlockHamRDeque.front(), herm(*CR)));
      }
   }

   std::deque<MatrixOperator>
   operator()(std::deque<MatrixOperator> const& XDeque) const
   {
      // Construct the "B"-matrices corresponding to the input "X"-matrices.
      std::deque<StateComponent> BDeque;
      auto NL = NullLeftDeque.begin();
      auto X = XDeque.begin();
      while (NL != NullLeftDeque.end())
      {
         BDeque.push_back(prod(*NL, *X));
         ++NL, ++X;
      }

      // Construct the "triangular" MPS unit cell.
      LinearWavefunction PsiTri;

      if (PsiLeft.size() == 1)
         PsiTri.push_back(BDeque.back());
      else
      {
         auto CL = PsiLinearLeft.begin();
         auto CR = PsiLinearRight.begin();
         auto B = BDeque.begin();
         SumBasis<VectorBasis> NewBasis0((*CL).Basis2(), (*B).Basis2());
         PsiTri.push_back(tensor_row_sum(*CL, *B, NewBasis0));
         ++CL, ++CR, ++B;
         for (int i = 1; i < PsiLeft.size()-1; ++i)
         {
            StateComponent Z = StateComponent((*CL).LocalBasis(), (*CR).Basis1(), (*CL).Basis2());
            SumBasis<VectorBasis> NewBasis1((*CL).Basis2(), (*B).Basis2());
            SumBasis<VectorBasis> NewBasis2((*CL).Basis1(), (*CR).Basis1());
            PsiTri.push_back(tensor_col_sum(tensor_row_sum(*CL, *B, NewBasis1), tensor_row_sum(Z, *CR, NewBasis1), NewBasis2));
            ++CL, ++CR, ++B;
         }
         SumBasis<VectorBasis> NewBasis3((*B).Basis1(), (*CR).Basis1());
         PsiTri.push_back(tensor_col_sum(*B, *CR, NewBasis3));
      }
      
      // Calcaulate the terms in the triangular E and F matrices where there is
      // one B-matrix on the top.
      StateComponent BlockHamLTri, BlockHamRTri;

      SolveSimpleMPO_Left2(BlockHamLTri, BlockHamL, PsiLinearLeft, PsiLinearRight, PsiTri,
                           PsiLeft.qshift(), HamMPO, RhoL, RhoL, ExpIK, GMRESTol, Verbose-1);

      SolveSimpleMPO_Right2(BlockHamRTri, BlockHamR, PsiLinearLeft, PsiLinearRight, PsiTri,
                            PsiRight.qshift(), HamMPO, RhoR, RhoR, ExpIK, GMRESTol, Verbose-1);

      // Shift the phases by one unit cell.
      BlockHamLTri *= ExpIK;
      BlockHamRTri *= ExpIK;

      // Calculate the contribution to HEff corresponding to where the
      // B-matrices are on the same site.
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
         ++B, ++NL, ++O, ++BHL, ++BHR;
      }

      // Calculate the contribution where the top B-matrix is in the left
      // semi-infinite part.
      StateComponent Tmp = BlockHamLTri;
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
         // TODO: We don't need to do this on the final step.
         Tmp = contract_from_left(*O, herm(*CR), Tmp, *CL) + contract_from_left(*O, herm(*B), *BHL, *CL);
         ++B, ++NL, ++CL, ++CR, ++O, ++BHL, ++BHR, ++R;
      }

      // Calculate the contribution where the top B-matrix is in the right
      // semi-infinite part.
      Tmp = BlockHamRTri;
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
         Tmp = contract_from_right(herm(*O), *CL, Tmp, herm(*CR)) + contract_from_right(herm(*O), *B, *BHR, herm(*CR));
      }

      return Result;
   }

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

   InfiniteWavefunctionLeft PsiLeft;
   InfiniteWavefunctionLeft PsiRight;
   BasicTriangularMPO HamMPO;
   double GMRESTol;
   int Verbose;
   std::complex<double> ExpIK;
   LinearWavefunction PsiLinearLeft, PsiLinearRight;
   StateComponent BlockHamL, BlockHamR;
   std::deque<StateComponent> BlockHamLDeque, BlockHamRDeque;
   std::deque<StateComponent> NullLeftDeque;
   MatrixOperator RhoL, RhoR;
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
      int NDisplay = 5;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("momentum,k", prog_opt::value(&k),
          FormatDefault("Excitation momentum (divided by pi)", k).c_str())
         ("ndisplay,n", prog_opt::value<int>(&NDisplay),
          FormatDefault("The number of lowest eigenvalues to display", NDisplay).c_str())
         ("maxiter", prog_opt::value<int>(&Iter),
          FormatDefault("Maximum number of Lanczos iterations", Iter).c_str())
         ("miniter", prog_opt::value<int>(&MinIter),
          FormatDefault("Minimum number of Lanczos iterations", MinIter).c_str())
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Error tolerance for the Lanczos eigensolver", Tol).c_str())
         ("gmrestol", prog_opt::value(&GMRESTol),
          FormatDefault("Error tolerance for the GMRES algorithm", GMRESTol).c_str())
         ("seed", prog_opt::value<unsigned long>(), "Random seed")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&InputFileLeft), "psi")
         ("ham", prog_opt::value(&HamStr), "ham")
         ("psi2", prog_opt::value(&InputFileRight), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("ham", 1);
      p.add("psi2", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("ham") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <hamiltonian> [psi-right]" << std::endl;
         std::cerr << desc << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX)
         : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      mp_pheap::InitializeTempPHeap();

      pvalue_ptr<MPWavefunction> InPsiLeft = pheap::ImportHeap(InputFileLeft);
      InfiniteWavefunctionLeft PsiLeft = InPsiLeft->get<InfiniteWavefunctionLeft>();

      InfiniteLattice Lattice;
      BasicTriangularMPO HamMPO, HamMPOLeft, HamMPORight;
      std::tie(HamMPO, Lattice) = ParseTriangularOperatorAndLattice(HamStr);

      // Rescale the momentum by the number of lattice unit cells in the unit cell of PsiLeft.
      k *= PsiLeft.size() / Lattice.GetUnitCell().size();

      HEff EffectiveHamiltonian;

      if (vm.count("psi2"))
      {
         pvalue_ptr<MPWavefunction> InPsiRight = pheap::ImportHeap(InputFileRight);
         InfiniteWavefunctionLeft PsiRight = InPsiRight->get<InfiniteWavefunctionLeft>();
         EffectiveHamiltonian = HEff(PsiLeft, PsiRight, HamMPO, k, GMRESTol, Verbose);
      }
      else
         EffectiveHamiltonian = HEff(PsiLeft, HamMPO, k, GMRESTol, Verbose);

      std::deque<MatrixOperator> XDeque = EffectiveHamiltonian.InitialGuess();

      LinearAlgebra::Vector<double> EValues = LanczosFull(XDeque, EffectiveHamiltonian, Iter, Tol, MinIter, Verbose);

      std::cout << "#n      En" << std::endl;

      auto E = EValues.begin();
      for (int i = 0; i < NDisplay && E != EValues.end(); ++i, ++E)
         std::cout << std::setw(4) << i << "    "
                   << std::setw(20) << *E << std::endl;
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
