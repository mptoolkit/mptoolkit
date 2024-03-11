// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// wavefunction/infinitewavefunctionleft.cpp
//
// Copyright (C) 2012-2024 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "infinitewavefunctionleft.h"
#include "infinitewavefunctionright.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "common/environment.h"

#include "common/statistics.h"

#include <fstream>

#include "wavefunction/operator_actions.h"
#include "pheap/pheapstream.h"
#include "interface/attributes.h"
#include "linearalgebra/arpack_wrapper.h"
#include "mps/packunpack.h"
#include "mp-algorithms/transfer.h"
#include "linearalgebra/takagi.h"
#include "linearalgebra/simdiag.h"

// Streaming versions:
// Note: the base class CanonicalWavefunctionBase has a separate version number.
//
// Version 1:
// No verioning information, this version needs to be deduced from the stream metadata version.
// Old style InfiniteWavefunction:
//      MatrixOperator        C_old
//      QuantumNumber         QShift
//      LinearWavefunctionOld PsiLinear
//      MatrixOperator        C_right
//      AttributeList         Attr
//
// Version 2:
//      CanonicalWavefunctionBase (base class)
//      QuantumNumber              QShift
//
// Version 3 adds:
//      double                     Amplitude
// Version 4:
//      remove 'Amplitude', replace by
//      double                     LogAmplitude

PStream::VersionTag
InfiniteWavefunctionLeft::VersionT(4);

double const Epsilon = std::numeric_limits<double>::epsilon();

extern double const ArnoldiTol = getenv_or_default("MP_ARNOLDI_TOL", Epsilon*10);

extern double const InverseTol = getenv_or_default("MP_INVERSE_TOL", 1E-7);

// the tol used in the orthogonalization can apparently be a bit smaller
extern double const OrthoTol = getenv_or_default("MP_ORTHO_TOL", 1E-8);

double const UnityEpsilon = 1E-14;

std::string InfiniteWavefunctionLeft::Type = "InfiniteWavefunctionLeft";

namespace
{

struct LeftMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiply(LinearWavefunction const& L_, QuantumNumber const& QShift_, int Verbose_ = 0)
      : L(L_), QShift(QShift_), Verbose(Verbose_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, QShift);
      int s = 0;
      for (LinearWavefunction::const_iterator I = L.begin(); I != L.end(); ++I)
      {
         if (Verbose > 0)
            std::cout << "site " << s << std::endl;
         r = operator_prod(herm(*I), r, *I);
         ++s;
      }
      return r;
   }

   LinearWavefunction const& L;
   QuantumNumber QShift;
   int Verbose;
};

struct RightMultiply
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiply(LinearWavefunction const& R_, QuantumNumber const& QShift_, int Verbose_ = 0)
      : R(R_), QShift(QShift_), Verbose(Verbose_) {}

   result_type operator()(argument_type const& x) const
   {
      result_type r = x;
      int s = R.size();
      LinearWavefunction::const_iterator I = R.end();
      while (I != R.begin())
      {
         --s;
         if (Verbose > 0)
            std::cout << "site " << s << std::endl;
         --I;
         r = operator_prod(*I, r, herm(*I));
      }
      return delta_shift(r, adjoint(QShift));
   }

   LinearWavefunction const& R;
   QuantumNumber QShift;
   int Verbose;
};

} // namespace

InfiniteWavefunctionLeft::InfiniteWavefunctionLeft(QuantumNumbers::QuantumNumber const& QShift_, double LogAmplitude_)
   : QShift(QShift_), LogAmplitude(LogAmplitude_)
{
}

InfiniteWavefunctionLeft::InfiniteWavefunctionLeft(InfiniteWavefunctionRight const& Psi)
{
   this->setBasis1(Psi.Basis1());
   this->setBasis2(Psi.Basis2());
   QShift = Psi.qshift();
   LogAmplitude = Psi.log_amplitude();

   // place-holder for the edge Lambda matrix (although in principle it doesn't change)
   this->push_back_lambda(Psi.lambda(0));

   auto PsiI = Psi.begin();
   auto LambdaI = Psi.lambda_begin();
   while (PsiI != Psi.end())
   {
      StateComponent A = (*LambdaI) * (*PsiI);
      MatrixOperator Vh;
      RealDiagonalOperator Lambda;
      std::tie(Lambda, Vh) = OrthogonalizeBasis2(A);
      A = A*Vh;
      MatrixOperator LambdaP = herm(Vh)*(Lambda*Vh);
      Lambda = ExtractRealDiagonal(LambdaP);
      // if (norm_frob(LambdaP-Lambda) > 1E-14)
      // {
      //    std::cerr << "get_right_canonical: warning: Lambda matrix isn't diagonal, difference = " << norm_frob(LambdaP-Lambda) << '\n';
      // }
      this->push_back(A);
      this->push_back_lambda(Lambda);
      ++PsiI;
      ++LambdaI;
   }
   this->set_lambda(0, delta_shift(this->lambda_r(), QShift));
   this->check_structure();
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::ConstructFromOrthogonal(LinearWavefunction Psi, QuantumNumbers::QuantumNumber const& QShift, RealDiagonalOperator const& Lambda, double LogAmplitude, int Verbose)
{
   InfiniteWavefunctionLeft Result(QShift, LogAmplitude);
   Result.InitializeFromLeftOrthogonal(Psi, Lambda, Verbose);
   return Result;
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::Construct(LinearWavefunction Psi, QuantumNumbers::QuantumNumber const& QShift, MatrixOperator GuessRho, double LogAmplitude, int Verbose)
{
   left_orthogonalize(Psi, QShift, ArnoldiTol, Verbose);
   auto Lambda = gauge_fix_left_orthogonal(Psi, QShift, GuessRho, ArnoldiTol, Verbose-1);
   return ConstructFromOrthogonal(Psi, QShift, Lambda, LogAmplitude, Verbose);
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::Construct(LinearWavefunction Psi,
                                    QuantumNumbers::QuantumNumber const& QShift_,
                                    double LogAmplitude,
                                    int Verbose)
{
   return InfiniteWavefunctionLeft::Construct(Psi, QShift_, MatrixOperator::make_identity(Psi.Basis2()), LogAmplitude, Verbose);
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::ConstructPreserveAmplitude(LinearWavefunction Psi, QuantumNumbers::QuantumNumber const& QShift, MatrixOperator GuessRho, double LogAmplitude, int Verbose)
{
   CHECK_EQUAL(GuessRho.Basis1(), Psi.Basis2());
   CHECK_EQUAL(GuessRho.Basis2(), Psi.Basis2());
   LogAmplitude += left_orthogonalize(Psi, QShift, ArnoldiTol, Verbose);
   auto Lambda = gauge_fix_left_orthogonal(Psi, QShift, GuessRho, ArnoldiTol, Verbose-1);
   return ConstructFromOrthogonal(Psi, QShift, Lambda, LogAmplitude, Verbose);
}

InfiniteWavefunctionLeft
InfiniteWavefunctionLeft::ConstructPreserveAmplitude(LinearWavefunction Psi, QuantumNumbers::QuantumNumber const& QShift_, double LogAmplitude, int Verbose)
{
   return InfiniteWavefunctionLeft::ConstructPreserveAmplitude(Psi, QShift_, MatrixOperator::make_identity(Psi.Basis2()), LogAmplitude, Verbose);
}

std::tuple<std::complex<double>, MatrixOperator>
GetPrincipalTransferEigenvectorLeft(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, double tol, int Verbose)
{
   std::vector<std::complex<double>> EValues;
   std::vector<MatrixOperator> EVec;

   std::tie(EValues, EVec) = get_left_transfer_eigenvectors(2, Psi, Psi, QShift, MatrixOperator::make_identity(Psi.Basis1()), tol, Verbose);

   // verify that the eigenvalues are non-degenerate
   if (EValues.size() > 1 && (std::abs(EValues[0] - EValues[1]) / (std::abs(EValues[0]) + std::abs(EValues[1])) <= UnityEpsilon))
   {
      std::cerr << "left_orthogonalize: error: largest eigenvalue of the transfer matrix is degenerate: "
                << EValues[0] << " and " << EValues[1] << '\n';
      std::abort(); //TODO: it would be better to throw an exception at this point
   }

   return std::tie(EValues[0], EVec[0]);
}

// helper function to return the minimum value of A[i]*sin(theta) + B[i]*cos(theta) with respect to i.
// We don't need the value of i, but just the minimum value itself.
double GetMinimumOfCombination(LinearAlgebra::Vector<double> const& A, LinearAlgebra::Vector<double> const& B, double Theta)
{
   double s = std::sin(Theta);
   double c = std::cos(Theta);
   double x = s*A[0] + c*B[0];
   for (int i = 1; i < A.size(); ++i)
   {
      double xx = s*A[i] + c*B[i];
      if (xx < x)
         x = xx;
   }
   return x;
}

#if 0
std::tuple<std::complex<double>, MatrixOperator, int>
GetPrincipalTransferEigenvectorLeftDegen(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, int Degen, double tol, int Verbose)
{
   std::vector<std::complex<double>> EValues;
   std::vector<MatrixOperator> EVec;

   std::tie(EValues, EVec) = get_left_transfer_eigenvectors(Degen+1, Psi, Psi, QShift, MatrixOperator::make_identity(Psi.Basis1()), tol, Verbose);

   // Keep increasing ExpectedDegeneracy and loop until we find a non-degenerate eigenvalue
   while (EValues.size() > 1 && (std::abs(EValues[0] - EValues.back()) / (std::abs(EValues[0]) + std::abs(EValues.back()) <= UnityEpsilon)))
   {
      ++Degen;
      if (Verbose > 0)
      {
         std::cerr << "GetPrincipalTransferEigenvectorLeftDegen: Eigenvalues are degenerate, increasing degeneracy to " << ExpectedDegeneracy << '\n';
      }
      std::tie(EValues, EVec) = get_left_transfer_eigenvectors(Degen+1, Psi, Psi, QShift, MatrixOperator::make_identity(Psi.Basis1()), tol, Verbose);
   }

   // Make sure the eigenvalues really are degenerate
   while (Degen > 1 && (std::abs(EValues[0] - EValues[Degen-1] / (std::abs(EValues[0]) + std::abs(EValues[Degen-1]) > UnityEpsilon))
      --Degen;

   // early return if the eigenvalues are non-degenerate
   if (Degen == 1)
   {
      return std::make_tuple(EValues[0], EVec[0], 1);
   }

   // Autonne–Takagi factorization of the degenerate eigenvectors
   LinearAlgebra::Matrix<std::complex<double>> A(Degen, Degen);
   for (int i = 0; i < Degen; ++i)
   {
      for (int j = 0; j < Degen; ++j)
      {
         A(i,j) = trace(EVec[i]*EVec[j]);
      }
   }
   LinearAlgebra::Matrix<std::complex<double>> U;
   LinearAlgebra::DiagonalMatrix<double> D;
   std::tie(U, D) = TakagiFactor(A);

   std::vector<MatrixOperator> X;
   for (int i = 0; i < Degen; ++i)
   {
      MatrixOperator ThisX = U(i,0) * EVec[0];
      for (int j = 1; j < Degen; ++j)
      {
         ThisX += U(i,j) * EVec[j];
      }
      X.push_back(std::move(ThisX));
   }

   LinearAlgebra::Matrix<std::complex<double>> U;
   std::vector<LinearAlgebra::Vector<double>> Y;
   std::tie(U, Y) = SimultaneousDiagonalizeHermitian(std::move(X));

   // We now have a diagonal set of matrices, and we want to find a positive linear combination.  This algorithm
   // only works for the Degen == 2 case.
   if (Degen > 2)
   {
      std::cerr << "GetPrincipalTransferEigenvectorLeftDegen: degeneracy is greater than 2, which is not yet supported.\n";
      std::abort();
   }
   double Theta = 0.0;
   double f = -1000;
   int N = size1(U);
   // try all possible values of theta
   for (int i = 0; i < N; ++i)
   {
      double TrialTheta = std::atan2(Y[0][i], Y[1][i]);
      double TrialF = GetMinimumOfCombination(Y[0], Y[1], TrialTheta);
      if (TrialF > f)
      {
         f = TrialF;
         Theta = TrialTheta;
      }
      // try the other quadrant
      if (TrialTheta <= 0)
         TrialTheta += math_const::pi;
      else
         TrialTheta -= math_const::pi;
      TrialF = GetMinimumOfCombination(Y[0], Y[1], TrialTheta);
      if (TrialF > f)
      {
         f = TrialF;
         Theta = TrialTheta;
      }
   }

   // Check that our final f is positive
   if (Verbose)
   {
      std::cerr << "Smallest element of the positive linear combination is " << f << '\n';
   }
   if (f < -1e-14)
   {
      std::cerr << "warning: smallest element of the positive linear combination is negative.\n";
   }

   MatrixOperator PositiveX = std::sin(Theta)*X[0] + std::cos(theta)*X[1];

   return std::tie(EValues[0], Positive, Degen);
}
#endif

// This is a helper function that does the hard work of left_orthogonalize.
double
left_orthogonalize_from_evector(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, std::complex<double> EVal, MatrixOperator X, int Verbose = 0)
{

   MatrixOperator U0;
   RealDiagonalOperator D0;
   std::tie(D0, U0) = DiagonalizeHermitian(X); // Now X = U^\dagger D U

   // Get the min and max elements of D0
   double m = 1.0;
   double M = 0.0;
   for (int i = 0; i < D0.Basis1().size(); ++i)
   {
      for (int j = 0; j < D0.Basis1().dim(i); ++j)
      {
         double x = D0(i,i)(j,j);
         m = std::min(m, x);
         M = std::max(M, x);
      }
   }
   double ConditionNum = M / m;
   if (Verbose > 0 || ConditionNum > 10.0)
   {
      std::cerr << "left_orthogonalize_from_evector: Left evector condition number is " << ConditionNum << '\n';
   }

   D0 = SqrtDiagonal(D0);

   MatrixOperator U = U0;
   RealDiagonalOperator D = D0;
   X = herm(U0) * (D*U);
   // Do a sequence of SVD's to orthogonalize Psi, and remove zero singular values
   for (auto& A : Psi)
   {
      A = X*A;
      std::tie(D, U) = TruncateBasis2(A);
      X = D*U;
   }
   // Incorporate the final unitaries.  This ensures that the final basis doesn't change.
   //Psi.set_front(herm(U0) * Psi.get_front()); // This would be the alternative; incorporate the unitary into the front instead
   Psi.set_back(Psi.get_back() * U);

   // Verify that the final eigenvector matrices are similar.  We need to scale by the eigenvalue.
   // Note that we compare D^2 here, since this is the relevant scale with respect to the numerical precision.
   MatrixOperator Dr = herm(U) * ((D*D) * U);
   Dr = delta_shift(Dr, QShift);
   MatrixOperator Dl = herm(U0) * ((D0*D0) * U0);
   Dl *= std::abs(EVal);

   double Eps = norm_frob(Dr-Dl);
   if (Verbose > 0)
   {
      std::cout << "Left transfer orthogonalization residual = " << Eps << '\n';
   }
   if (Eps >= 1E-12*Dr.Basis1().total_dimension())
   {
      std::cerr << "left_orthogonalize_from_evector: WARNING: left transfer orthogonalization residual is large: "
         << Eps << '\n';
   }

   // scale by 0.5 here to give the scaling for the A-matrices, since the overlap gives |amplitude|^2
   return 0.5*std::log(std::abs(EVal));
}

double
left_orthogonalize(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, double tol, int Verbose)
{
   // Strategy:
   // 1. Calculate the two largest left transfer matrix eigenvalues.
   // 2. Check and see if they are degenerate.  If they are, use the degeneracy algorithm.
   // 3. (non-degenerate case): Scale the eigenmatrix X to be hermitian
   // 4. Diagonalize it: X = U D U^\dagger
   // 5. Obtain the decompositon L = U \sqrt(D) so that L L^\dagger = X
   // 6. Do a sequence of SVD's from left to right
   // 7. Check that the final D matrix matches the initial

   std::complex<double> EVal;
   MatrixOperator X;
   // std::tie(EVal, X) = GetPrincipalTransferEigenvectorLeft(Psi, QShift, tol, Verbose-2);
   // std::tie(EVal, X) = get_left_transfer_eigenvector_hermitian(Psi, QShift, tol, Verbose-2);
   std::tie(EVal, X) = get_left_transfer_eigenvector(Psi, Psi, QShift, tol, Verbose-2);

   MatrixOperator XH = MatrixOperator::make_identity(X.Basis1()) * herm(X);
   X = 0.5*(X + XH);

   std::tie(EVal, X) = get_left_transfer_eigenvector(Psi, Psi, QShift, X, tol, Verbose-2);

   if (Verbose >= 2 && std::abs(EVal-1.0) > UnityEpsilon)
   {
      std::cerr << "left_orthogonalize: WARNING: left transfer eigenvalue differs from 1.0: "
                << EVal << " 1-evalue=" << (EVal-1.0) << '\n';
   }

   // scale X so that it is Hermitian, and positive
   auto x = trace(X*X);   // TODO: this isn't as efficient as it should be
   X *= std::pow(x, -0.5);
   if (trace(X).real() < 0)
      X = -X;

   // See how close X is to being Hermitian.  If it is not Hermitian, then we are probably
   // close to having degenerate eigenvalues.  We can try rerunning the eigensolver to get a
   // better eigenvector.
   MatrixOperator Xh = MatrixOperator::make_identity(X.Basis1()) * herm(X);
   double e = norm_frob(X-Xh);
   X = 0.5*(X+Xh);
   int Tries = 0;
   while (e > 1E-14 && ++Tries < 10) // some arbitrary number of tries.  In practice nothing seems to change after 2
   {
      std::cerr << "left_orthogonalize: warning: left transfer eigenmatrix is not hermitian, norm of X-X^dagger = " << e << " restarting orthogonalization\n";

      std::tie(EVal, X) = get_left_transfer_eigenvector(Psi, Psi, QShift, X);
      x = trace(X*X);   // TODO: this isn't as efficient as it should be
      X *= std::pow(x, -0.5);
      if (trace(X).real() < 0)
      X = -X;
      Xh = MatrixOperator::make_identity(X.Basis1()) * herm(X);
      e = norm_frob(X-Xh);
      X = 0.5*(X+Xh);
   }
   if (Tries == 10)
   {
      std::cerr << "left_orthogonalize: warning: left transfer eigenmatrix never converged to a Hermitian matrix!\n";
   }

   return left_orthogonalize_from_evector(Psi, QShift, EVal, X, Verbose);
}

#if 0
std::tuple<double, int>
left_orthogonalize_degen(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, int ExpectedDegen, double tol, int Verbose)
{
   std::complex<double> EVal;
   MatrixOperator X;
   int Degen,
   std::tie(EVal, X, Degen) = GetPrincipalTransferEigenvectorLeftDegen(Psi, QShift, ExpectedDegen, tol, Verbose-2);
   return std::make_tuple(left_orthogonalize_from_evector(Psi, QShift, EVal, X), Degen);
}

std::vector<LinearWavefunction>
left_branch_degen(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, int Degeneracy, double tol, int Verbose );
{
   std::vector<std::complex<double>> EValues;
   std::vector<MatrixOperator> EVec;

   std::tie(EValues, EVec) = get_left_transfer_eigenvectors(Degen, Psi, Psi, QShift, MatrixOperator::make_identity(Psi.Basis1()), tol, Verbose);

   // Autonne–Takagi factorization of the degenerate eigenvectors to make a hermitian basis
   LinearAlgebra::Matrix<std::complex<double>> A(Degen, Degen);
   for (int i = 0; i < Degen; ++i)
   {
      for (int j = 0; j < Degen; ++j)
      {
         A(i,j) = trace(EVec[i]*EVec[j]);
      }
   }
   LinearAlgebra::Matrix<std::complex<double>> U;
   LinearAlgebra::DiagonalMatrix<double> D;
   std::tie(U, D) = TakagiFactor(A);

   std::vector<MatrixOperator> X;
   for (int i = 0; i < Degen; ++i)
   {
      MatrixOperator ThisX = U(i,0) * EVec[0];
      for (int j = 1; j < Degen; ++j)
      {
         ThisX += U(i,j) * EVec[j];
      }
      X.push_back(std::move(ThisX));
   }

   // Since the MPS is left orthogonal, we know that one eigenvector is the identity.  So we can orthogonalize the
   // X matrices against the identity and there will be N-1 linearly independent results.  In the case where
   // N=2, this is easlier.
}
#endif

RealDiagonalOperator
gauge_fix_left_orthogonal(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, MatrixOperator GuessRho, double tol, int Verbose)
{
   // 8. Calculate the right transfer matrix eigenvector (density operator) \rho
   // 9. Diagonalize \rho and obtain U_L \Lambda_L^2 U_L^\dagger
   // 10. Gauge fix the left side with U_L^\dagger
   // 11. We now have a left-orthogonal MPS, so we can ConstructFromLeftOrthogonal

   CHECK_EQUAL(GuessRho.Basis1(), Psi.Basis2());
   CHECK_EQUAL(GuessRho.Basis2(), Psi.Basis2());

   std::complex<double> EValue;
   MatrixOperator Y;
   // Ensure the initial guess is hermitian
   MatrixOperator rH = MatrixOperator::make_identity(GuessRho.Basis1()) * herm(GuessRho);
   GuessRho = 0.5 * (GuessRho + rH);
   // std::tie(EValue, Y) = get_right_transfer_eigenvector_hermitian(Psi, QShift, std::move(GuessRho), tol, Verbose);
   std::tie(EValue, Y) = get_right_transfer_eigenvector(Psi, Psi, QShift, std::move(GuessRho), tol, Verbose);
   if (std::abs(EValue-1.0) > UnityEpsilon)
   {
      std::cerr << "gauge_fix_left_orthogonal: WARNING: right transfer eigenvalue differs from 1.0: "
      << EValue << " evalue-1=" << (EValue-1.0) << '\n';
   }
   // scale Y so it is Hermitian - as an eigenvector it may be arbitrarily rotated by a complex number
   auto y = trace(Y*Y);  // TODO: fixme
   Y *= std::pow(y, -0.5);
   // make the trace of Y = 1
   Y *= 1.0 / trace(Y);
   // if (trace(Y).real() < 0.0)
   //    Y *= -1.0;

   rH = MatrixOperator::make_identity(Y.Basis1()) * herm(Y);
   double e = norm_frob(Y-rH);
   Y = 0.5*(Y+rH);
   int Tries = 0;
   while (e > 1E-14 && ++Tries < 10)
   {
      std::cerr << "gauge_fix_left_orthogonal: warning: right transfer eigenmatrix is not hermitian, norm of Y-Y^dagger = " << e << " restarting orthogonalization\n";
      std::tie(EValue, Y) = get_right_transfer_eigenvector(Psi, Psi, QShift, std::move(Y), tol, Verbose);
      // scale Y so it is Hermitian - as an eigenvector it may be arbitrarily rotated by a complex number
      auto y = trace(Y*Y);  // TODO: fixme
      Y *= std::pow(y, -0.5);
      // make the trace of Y = 1
      Y *= 1.0 / trace(Y);

      rH = MatrixOperator::make_identity(Y.Basis1()) * herm(Y);
      e = norm_frob(Y-rH);
      Y = 0.5*(Y+rH);
   }
   if (Tries == 10)
   {
      std::cerr << "gauge_fix_left_orthogonal: warning: right transfer eigenmatrix never converged to a Hermitian matrix!\n";
   }

   MatrixOperator U;
   RealDiagonalOperator D;
   std::tie(D, U) = DiagonalizeHermitian(std::move(Y)); // Now Y = U^\dagger D U
   // FIXME: we probably should remove zero singular values from D here, if there are any

   Psi.set_front(delta_shift(U, QShift) * Psi.get_front());
   Psi.set_back(Psi.get_back() * herm(U));

   return SqrtDiagonal(std::move(D));
}

RealDiagonalOperator
gauge_fix_left_orthogonal(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift, double tol, int Verbose)
{
   return gauge_fix_left_orthogonal(Psi, QShift, MatrixOperator::make_identity(Psi.Basis1()), tol, Verbose);
}

void
InfiniteWavefunctionLeft::InitializeFromLeftOrthogonal(LinearWavefunction Psi, RealDiagonalOperator Lambda, int Verbose)
{
   if (Verbose > 1)
   {
      std::cout << "Constructing Lambda matrices..." << std::endl;
   }

   MatrixOperator U1 = MatrixOperator::make_identity(Lambda.Basis1());

   RealDiagonalOperator Lambda0 = delta_shift(Lambda, this->QShift);

   std::deque<RealDiagonalOperator> LambdaMatrices;
   LambdaMatrices.push_front(Lambda);

   LinearWavefunction::iterator I = Psi.end();
   while (I != Psi.begin())
   {
      --I;

      StateComponent A = *I;

      A = A * (U1 * Lambda);
      MatrixOperator U0;
      std::tie(U0, Lambda) = OrthogonalizeBasis1(A);

      // We incorporate the left-most U into the right basis, which preserves the basis at the unit cell boundary
      if (I != Psi.begin())
         *I = triple_prod(herm(U0), *I, U1);
      else
         *I = (*I) * U1;
      LambdaMatrices.push_front(Lambda);
      U1 = U0;
   }

   // Now we have a choice as to which basis we use at the edge of the unit cell.  We can either use
   // Lambda, which means we need to incorporate U1 at the right edge, or we could use Lambda0, and
   // incorporate U1 at the left edge.  We will use Lambda0, since that preserves the existing basis.
   LambdaMatrices.front() = Lambda0;

   // check that Lambda0 and Lambda are the same (up to a delta_shift, and epsilon).
   MatrixOperator L0 = Lambda0*Lambda0;
   MatrixOperator L = U1 * (Lambda*Lambda) * herm(U1);
   double Eps = norm_frob(L0 - L);
   if (Verbose > 0)
   {
      std::cout << "Right transfer orthogonalization residual = " << Eps << '\n';
   }
   if (Eps >= Lambda.Basis1().total_dimension() * std::sqrt(Psi.size()) * 1E-14)
   {
      std::cerr << "InfiniteWavefunctionLeft: WARNING: right transfer orthogonalization residual is large: "
         << Eps << '\n';
      TRACE(Lambda0)(Lambda);
   }

   this->set_A_matrices_from_handle(Psi.base_begin(), Psi.base_end());
   this->set_lambda_matrices(LambdaMatrices.begin(), LambdaMatrices.end());

   this->setBasis1(Psi.Basis1());
   this->setBasis2(Psi.Basis2());

   this->check_structure();
}

void read_version(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi, int Version)
{
   if (Version == 1)
   {
      PANIC("This wavefunction is too old.");
      // ConstructFromOrthogonal now takes a RealDiagonalOperator.
      // So we just construct with C_right*C_right as an initial guess for the density matrix.
      MatrixOperator C_old;
      QuantumNumbers::QuantumNumber QShift;
      LinearWavefunction PsiLinear;
      MatrixOperator C_right;
      AttributeList Attr;
      AttributeList AttrDummy;

      in >> C_old;
      in >> QShift;
      in >> PsiLinear;
      in >> AttrDummy;
      in >> C_right;
      in >> Attr;

      Psi = InfiniteWavefunctionLeft::Construct(PsiLinear, QShift, C_right*C_right);
   }
   else
   {
      int BaseVersion = Psi.CanonicalWavefunctionBase::ReadStream(in);
      in >> Psi.QShift;
      if (BaseVersion < 3)
         Psi.set_lambda(0, delta_shift(Psi.lambda_r(), Psi.qshift()));

      if (Version == 3)
      {
         // new in version 3, but replaced by LogAmplitude in version 4
         double Amplitude;
         in >> Amplitude;
         Psi.LogAmplitude = std::log(Amplitude);
      }
      else if (Version >= 4)
      {
         in >> Psi.LogAmplitude;
      }
      else
      {
         // prior to version 3, the amplitude is assumed to be always 1
         Psi.LogAmplitude = 0.0;
      }

      if (Version > 4)
      {
         PANIC("This program is too old to read this wavefunction, expected Version <= 4")(Version);
      }
   }

   Psi.debug_check_structure();
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteWavefunctionLeft& Psi)
{
   // We formerly didn't have a version number for InfiniteWavefunction, so we have a hack,
   // and read the metadata version of the filesystem.  This is 1 if the file is old and
   // doesn't contain an InfiniteWavefunction version number.
   int Version;

   PHeapFileSystem::ipheapstream* S = dynamic_cast<PHeapFileSystem::ipheapstream*>(&in);
   if (S && S->version() == 1)
   {
      // old file, need to set the version number manually
      Version = 1;
   }
   else
   {
      // new file, read the version number
      Version = in.read<int>();
   }

   PStream::VersionSentry Sentry(in, InfiniteWavefunctionLeft::VersionT, Version);
   read_version(in, Psi, Version);

   return in;
}

PStream::opstream& operator<<(PStream::opstream& out, InfiniteWavefunctionLeft const& Psi)
{
   out << InfiniteWavefunctionLeft::VersionT.default_version();

   Psi.CanonicalWavefunctionBase::WriteStream(out);

   out << Psi.QShift;
   out << Psi.LogAmplitude;

   return out;
}

std::pair<LinearWavefunction, RealDiagonalOperator>
get_left_canonical(InfiniteWavefunctionLeft const& Psi)
{
   return std::make_pair(LinearWavefunction(Psi.Basis1().GetSymmetryList(),
                                            Psi.base_begin(), Psi.base_end()), Psi.lambda_r());
}

InfiniteWavefunctionLeft
coarse_grain(InfiniteWavefunctionLeft const& Psi, int N, int Verbose)
{
   // TODO: this could be improved, because we don't need to reconstruct all of the Lambda matrices.
   LinearWavefunction PsiL;
   RealDiagonalOperator Lambda;
   std::tie(PsiL, Lambda) = get_left_canonical(Psi);
   return InfiniteWavefunctionLeft::ConstructFromOrthogonal(coarse_grain(PsiL, N), Psi.qshift(), Lambda, Psi.log_amplitude()*N, Verbose);
}

std::tuple<RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionLeft const& Psi)
{
   // this could be parallelised to do all sites at once
   LinearWavefunction Result;
   RealDiagonalOperator D;
   auto PsiI = Psi.begin();
   auto LambdaI = Psi.lambda_begin();
   ++LambdaI; // skip to the first bond within the unit cell
   while (PsiI != Psi.end())
   {
      StateComponent A = (*PsiI) * (*LambdaI);
      MatrixOperator U;
      RealDiagonalOperator Lambda;
      std::tie(U, Lambda) = OrthogonalizeBasis1(A);
      A = U*A;
      MatrixOperator LambdaP = U*Lambda*herm(U);
      Lambda = ExtractRealDiagonal(LambdaP);
      // if (norm_frob(LambdaP-Lambda) > 1E-14)
      // {
      //    std::cerr << "get_right_canonical: warning: Lambda matrix isn't diagonal, difference = " << norm_frob(LambdaP-Lambda) << '\n';
      // }
      Result.push_back(A);
      if (PsiI == Psi.begin())
         D = Lambda;
      ++PsiI;
      ++LambdaI;
   }
   return std::make_tuple(D, Result);
}

void InfiniteWavefunctionLeft::scale(std::complex<double> x)
{
   double Amplitude = std::abs(x);
   x = x / Amplitude;  // x is now the phase
   LogAmplitude += std::log(Amplitude);
   // Because the components have value semantics we can't directly load() the pvalue_ptr
   // and then mutate it, because pvalue_handle's are immutable - we would end up referring
   // do a different object.
   pvalue_ptr<StateComponent> s = this->base_begin_()->load();
   *s.mutate() *= x;
   *this->base_begin_() = s;
}

void InfiniteWavefunctionLeft::scale_log(std::complex<double> x)
{
   LogAmplitude += x.real();
   // Because the components have value semantics we can't directly load() the pvalue_ptr
   // and then mutate it, because pvalue_handle's are immutable - we would end up referring
   // do a different object.
   if (x.imag() != 0.0)
   {
      pvalue_ptr<StateComponent> s = this->base_begin_()->load();
      *s.mutate() *= std::exp(std::complex<double>(0,x.imag()));
      *this->base_begin_() = s;
   }
}

void InfiniteWavefunctionLeft::normalize()
{
   LogAmplitude = 0;
}

InfiniteWavefunctionLeft& operator*=(InfiniteWavefunctionLeft& psi, std::complex<double> x)
{
   psi.scale(x);
   return psi;
}

void
InfiniteWavefunctionLeft::rotate_left(int Count)
{
   // Rotation is fairly straightforward, we just rotate the vectors around
   if (Count < 0)
   {
      this->rotate_right(-Count);
      return;
   }

   Count = Count % this->size();
   if (Count == 0)
      return;

   // the first Count elements are going to get shifted to the right hand side, so we need to
   // delta_shift them
   for (mps_iterator I = this->begin_(); I != this->begin_()+Count; ++I)
   {
      I->delta_shift(adjoint(this->qshift()));
   }
   // now do the actual rotation
   std::rotate(this->base_begin_(), this->base_begin_()+Count, this->base_end_());

   // for the Lambda matrices, start by removing the double-counted boundary lambda
   this->pop_back_lambda();
   // and delta-shift
   for (lambda_iterator I = this->lambda_begin_(); I != this->lambda_begin_()+Count; ++I)
   {
      I->delta_shift(adjoint(this->qshift()));
   }
   // and rotate
   std::rotate(this->lambda_base_begin_(), this->lambda_base_begin_()+Count, this->lambda_base_end_());
   // and put back the boundary lambda
   this->push_back_lambda(delta_shift(this->lambda_l(), adjoint(this->qshift())));

   // set the left and right basis
   this->setBasis1(lambda_l().Basis1());
   this->setBasis2(lambda_r().Basis2());

   this->debug_check_structure();
}

void
InfiniteWavefunctionLeft::rotate_right(int Count)
{
   if (Count < 0)
   {
      this->rotate_left(-Count);
      return;
   }

   Count = Count % this->size();
   if (Count == 0)
      return;

   this->rotate_left(this->size() - Count);
}

void
InfiniteWavefunctionLeft::check_structure() const
{
   this->CanonicalWavefunctionBase::check_structure();

   CHECK_EQUAL(this->Basis1(), delta_shift(this->Basis2(), this->qshift()));
}

void inplace_reflect(InfiniteWavefunctionLeft& Psi)
{
   // Although we construct the reflection into a new wavefunction,
   // we do gain an advantage for the inplace operation because we remove
   // elements from Psi as we go, so memory and disk use is ~ constant
   InfiniteWavefunctionLeft Result;
   Result.QShift = Psi.qshift();

   int Size = Psi.size();

   RealDiagonalOperator D = Psi.lambda_r();
   MatrixOperator U = MatrixOperator::make_identity(D.Basis1());

   // left-most lambda matrix, this is a place-holder
   Result.push_back_lambda(flip_conj(D));
   RealDiagonalOperator DSave = D;

   InfiniteWavefunctionLeft::const_base_mps_iterator I = Psi.base_end();
   while (I != Psi.base_begin())
   {
      --I;

      StateComponent A = prod(*I->lock(), U*D);
      MatrixOperator M = ExpandBasis1(A);

      MatrixOperator Vh;
      SingularValueDecomposition(M, U, D, Vh);
      A = prod(Vh, A);

      Result.push_back(reflect(A));
      Result.push_back_lambda(flip_conj(D));

      Psi.pop_back();
      Psi.pop_back_lambda();
   }

   // We have the final U matrix to deal with.  U.Basis1() is our original basis,
   // U.Basis2() has some gauge transform.  We want to write the final wavefunction in the
   // (flip conjugate) of the original basis.  If the wavefunction is perfectly orthogonal then
   // we should satisfy DSave = U*D*herm(U) here.  We want to use the DSave basis, not the D basis.
   // So we need to form
   // U*D (from the final SVD) = U*D*herm(U)*U = DSave * U
   // and merge the U into the final A-matrix.

#if 1
   //   Result.set_lambda(0, delta_shift(flip_conj(DSave), Psi.qshift()));
   //   Result.set_lambda(Size, flip_conj(DSave));
   Result.set_lambda(0, flip_conj(DSave));
   Result.set_lambda(Size, delta_shift(flip_conj(DSave), adjoint(Psi.qshift())));

   Result.set(Size-1, prod(Result[Size-1], herm(flip_conj(U))));

   Result.setBasis1(Result[0].Basis1());
   Result.setBasis2(Result[Size-1].Basis2());
   Result.LogAmplitude = Psi.LogAmplitude;
#else
   // old code that used the D basis (and hence introduces a gauge transformation)
   Result.set_lambda(0, delta_shift(flip_conj(D), Psi.qshift()));
   Result.set(0, prod(herm(delta_shift(flip_conj(U), Psi.qshift())), Result[0]));

   Result.setBasis1(Result[0].Basis1());
   Result.setBasis2(adjoint(D.Basis2()));
#endif

   Psi = Result;

   Psi.check_structure();
}

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionLeft& Psi)
{
   for (InfiniteWavefunctionLeft::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = conj(*I);
   }
}

void inplace_qshift(InfiniteWavefunctionLeft& Psi, QuantumNumbers::QuantumNumber const& Shift)
{
   Psi.setBasis1(delta_shift(Psi.Basis1(), Shift));
   Psi.setBasis2(delta_shift(Psi.Basis2(), Shift));

   for (InfiniteWavefunctionLeft::mps_iterator I = Psi.begin_(); I != Psi.end_(); ++I)
   {
      *I = delta_shift(*I, Shift);
   }

   for (InfiniteWavefunctionLeft::lambda_iterator I = Psi.lambda_begin_(); I != Psi.lambda_end_(); ++I)
   {
      *I = delta_shift(*I, Shift);
   }

   Psi.check_structure();
}

InfiniteWavefunctionLeft repeat(InfiniteWavefunctionLeft const& Psi, int Count)
{
   CHECK(Count >= 0);

   if (Count == 0 || Psi.empty())
      return InfiniteWavefunctionLeft();

   if (Count == 1)
      return Psi;

   // Claculate the final QShift, we need this for later
   QuantumNumber FinalQShift(Psi.GetSymmetryList());
   for (int i = 0; i < Count; ++i)
      FinalQShift = delta_shift(FinalQShift, Psi.qshift());

   // For the lambda matrices, we don't want to add the boundaries twice
   InfiniteWavefunctionLeft::const_lambda_iterator LambdaE = Psi.lambda_end();
   --LambdaE;

   // we need to handle the delta shift
   InfiniteWavefunctionLeft Result;
   bool First = true;
   for ( ; Count > 1; --Count)
   {
      QuantumNumber q(Psi.GetSymmetryList());
      for (int i = 1; i < Count; ++i)
         q = delta_shift(q, Psi.qshift());

      for (InfiniteWavefunctionLeft::const_mps_iterator I = Psi.begin(); I != Psi.end(); ++I)
      {
         // The very first time through the loop, set the Basis1.
         // We are guaranteed to go through this loop at least once, because
         // we return early if Count < 2
         if (First)
         {
            StateComponent A = delta_shift(*I, q);
            Result.setBasis1(A.Basis1());
            Result.push_back(A);
            Result.QShift = q;
            First = false;
         }
         else
         {
            Result.push_back(delta_shift(*I, q));
         }
      }

      for (InfiniteWavefunctionLeft::const_lambda_iterator I = Psi.lambda_begin(); I != LambdaE; ++I)
      {
         Result.push_back_lambda(delta_shift(*I, q));
      }

   }

   // The last repetition can make a simple copy
   for (InfiniteWavefunctionLeft::const_base_mps_iterator I = Psi.base_begin(); I != Psi.base_end(); ++I)
   {
      Result.push_back(*I);
   }

   // Here, we do want to go right to the end
   for (InfiniteWavefunctionLeft::const_base_lambda_iterator I = Psi.lambda_base_begin();
        I != Psi.lambda_base_end(); ++I)
   {
      Result.push_back_lambda(*I);
   }


   Result.setBasis2(Psi.Basis2());
   Result.QShift = FinalQShift;
   Result.LogAmplitude = Count * Psi.log_amplitude();

   return Result;
}

std::tuple<std::vector<std::complex<double>>, int>
overlap(InfiniteWavefunctionLeft const& x, ProductMPO const& StringOp,
	       InfiniteWavefunctionLeft const& y, int NumEigen,
	       QuantumNumbers::QuantumNumber const& Sector, bool UseAmplitude, double Tol, int Verbose)
{
   int Length = statistics::lcm(x.size(), y.size(), StringOp.size());

   ProductMPO Str = StringOp * ProductMPO::make_identity(StringOp.LocalBasis2List(), Sector);

   LinearWavefunction xPsi = get_left_canonical(x).first;
   LinearWavefunction yPsi = get_left_canonical(y).first;

   ApplyToPackedStateComponent<LeftMultiplyOperator> Func(LeftMultiplyOperator(xPsi, x.qshift(), Str, yPsi, y.qshift(), Length),
							  Str.Basis1(), x.Basis1(), y.Basis1());
   int n = Func.pack_size();
   int ncv = 0;

   LinearAlgebra::Vector<std::complex<double>> Eigen =
      LinearAlgebra::DiagonalizeARPACK(Func, n, NumEigen, nullptr, Tol, nullptr, ncv, true, Verbose);

	// Scale the eigenvalues by the amplitude
	if (UseAmplitude)
	{
		double LogAmplitude = x.log_amplitude() * (Length/x.size()) + y.log_amplitude() * (Length/y.size());
		for (auto& e : Eigen)
			e *= std::exp(LogAmplitude);
	}

   return std::make_tuple(std::vector<std::complex<double>>(Eigen.begin(), Eigen.end()), Length);
}


InfiniteWavefunctionLeft
reflect(InfiniteWavefunctionRight const& Psi)
{
   PANIC("not implemented");
}

void
InfiniteWavefunctionLeft::SetDefaultAttributes(AttributeList& A) const
{
   A["WavefunctionType"] = "Infinite";
   A["UnitCellSize"] = this->size();
   A["QShift"] = this->qshift();
}

// inject_left for a BasicFiniteMPO.  This can have support on multiple wavefunction unit cells
MatrixOperator
inject_left(MatrixOperator const& m,
            InfiniteWavefunctionLeft const& Psi1,
            BasicFiniteMPO const& Op,
            InfiniteWavefunctionLeft const& Psi2)
{
   CHECK_EQUAL(Psi1.size(), Psi2.size());
   CHECK_EQUAL(Psi1.qshift(), Psi2.qshift());
   DEBUG_CHECK_EQUAL(m.Basis1(), Psi1.Basis1());
   DEBUG_CHECK_EQUAL(m.Basis2(), Psi2.Basis1());
   if (Op.is_null())
      return MatrixOperator();

   // we currently only support simple irreducible operators
   CHECK_EQUAL(Op.Basis1().size(), 1);
   CHECK_EQUAL(Op.Basis2().size(), 1);
   CHECK_EQUAL(Op.Basis1()[0], m.TransformsAs());
   MatrixOperator Result = m;
   StateComponent E(Op.Basis1(), m.Basis1(), m.Basis2());
   E[0] = m;
   E.debug_check_structure();
   InfiniteWavefunctionLeft::const_mps_iterator I1 = Psi1.begin();
   InfiniteWavefunctionLeft::const_mps_iterator I2 = Psi2.begin();
   BasicFiniteMPO::const_iterator OpIter = Op.begin();
   while (OpIter != Op.end())
   {
      if (I1 == Psi1.end())
      {
         I1 = Psi1.begin();
         I2 = Psi2.begin();
         E = delta_shift(E, Psi1.qshift());
      }
      E = contract_from_left(*OpIter, herm(*I1), E, *I2);
      ++I1; ++I2; ++OpIter;
   }
   return delta_shift(E[0], Psi1.qshift());
}

std::complex<double>
expectation(InfiniteWavefunctionLeft const& Psi, BasicFiniteMPO const& Op)
{
#if 0
   MatrixOperator X = MatrixOperator::make_identity(Psi.Basis1());
   X = inject_left(X, Psi, Op, Psi);

   MatrixOperator Rho = Psi.lambda_r();
   Rho = Rho*Rho;

   return inner_prod(delta_shift(Rho, Psi.qshift()), X);
#else
   // find the first non-trivial site of the MPO
   BasicFiniteMPO::const_iterator i = Op.begin();
   int n = 0;
   OperatorClassification c;
   c.Identity_ = true;
   if (i != Op.end())
      c = classify(*i);
   while (i != Op.end() && c.is_identity())
   {
      ++i; ++n;
      if (i != Op.end())
         c = classify(*i);
   }
   // find the last non-trivial site of the MPO
   BasicFiniteMPO::const_iterator j = Op.end();
   int m = Op.size();
   c = OperatorClassification();
   c.Identity_ = true;
   while (j != i && c.is_identity())
   {
      --j; --m;
      c = classify(*j);
   }
   if (i == j && c.is_identity())
      return 1.0;
   // we could optimize for the case i==j and c.is_identity()
   ++j; ++m;
   // get an iterator into the wavefunction
   InfiniteWavefunctionLeft::const_mps_iterator I = Psi.begin();
   for (int a = 0; a < n; ++a)
   {
      ++I;
      if (I == Psi.end())
         I = Psi.begin();
   }
   //TRACE(m-n);
   MatrixOperator X = MatrixOperator::make_identity(I->Basis1());
   CHECK_EQUAL(i->Basis1().size(), 1);
   CHECK(is_scalar(i->Basis1()[0]));
   StateComponent E(i->Basis1(), X.Basis1(), X.Basis2());
   E[0] = X;
   while (n < m)
   {
      if (I == Psi.end())
      {
         I = Psi.begin();
         E = delta_shift(E, Psi.qshift());
      }
      E = contract_from_left(*i, herm(*I), E, *I);
      ++I;
      ++i;
      ++n;
   }
   m = m % Psi.size();
   if (m == 0)
      m = Psi.size();    // if we're at the edge of the unit cell, use lambda_r
   MatrixOperator Rho = Psi.lambda(m);
   Rho = Rho*Rho;

   return inner_prod(Rho, E[0]);
#endif
}
