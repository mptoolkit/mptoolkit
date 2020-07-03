// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/tebd.cpp
//
// Copyright (C) 2020 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mp-algorithms/tebd.h"
#include "mps/density.h"

// returns the number of digits of precision used in the decimal number f
int Digits(std::string const& f)
{
   std::string::const_iterator dot = std::find(f.begin(), f.end(), '.');
   if (dot == f.end())
      return 0;
   ++dot;
   std::string::const_iterator x = dot;
   while (x != f.end() && isdigit(*x))
      ++x;
   if (x != f.end() && (*x == 'e' || *x == 'E'))
   {
      // exponential notation
      ++x;
      int Exp = std::stoi(std::string(x, f.end()));
      return std::max(0, int((x-dot) - Exp));
   }
   return x-dot;
}

std::string FormatDigits(double x, int Digits)
{
   std::ostringstream str;
   str << std::fixed << std::setprecision(Digits) << std::setfill('0') << x;
   return str.str();
}

// Symmetric decompositions
// Even slice A, odd slice B.
// We always have either a.size() == b.size() OR a.size() == b.size()+1.
// If a.size() == b.size()+1 = n+1, then the pattern is A_1 B_1 A_2 B_2 ... A_N B_N A_{N+1} B_N A_N .... B_2 A_2 B_1 A_1
// if b.size() == a.size() = n, then the pattern is A_1 B_1 A_2 B_2 ... A_N B_N A_{N-1} B_{N-2} A_{N-2} .... B_1 A_1
// We don't need to include the final term of each array on construction since it is implicit:
// If a.size() == b.size()+1 then add final term to a of 1.0 - 2*sum(a) and final term to b of 0.5 - sum(b)
// If a.size() == b.size() then add final term to a of 0.5 - sum(a), and final term to b of 1.0 - 2*sum(b)

LTSDecomposition
SymmetricDecomposition(int Order, std::string Description,
                       std::vector<double> a,
                       std::vector<double> b)
{
   if (a.size() == b.size())
   {
      a.push_back(0.5 - std::accumulate(a.begin(), a.end(), 0.0));
      b.push_back(1.0 - 2.0 * std::accumulate(b.begin(), b.end(), 0.0));

      int asz = a.size();
      for (int i = asz-1; i >=0; --i)
         a.push_back(a[i]);

      int bsz = b.size();
      for (int i = bsz-2; i >= 0; --i)
         b.push_back(b[i]);
   }
   else if (a.size() == b.size()+1)
   {
      a.push_back(1.0 - 2.0 * std::accumulate(a.begin(), a.end(), 0.0));
      b.push_back(0.5 - std::accumulate(b.begin(), b.end(), 0.0));

      int asz = a.size();
      for (int i = asz-2; i >=0; --i)
         a.push_back(a[i]);

      int bsz = b.size();
      for (int i = bsz-1; i >= 0; --i)
         b.push_back(b[i]);
   }
   else
   {
      PANIC("Invalid SymmetricDecomposition");
   }

   return LTSDecomposition(Order, Description, a, b);
}

// The LeapfrogDecomposition assumes that the number of terms (m-1)/2 is odd, otherwise we need to repeat
// the final term again. We only need to supply n-1 weights, since the final weight is always
// 1 - 2 * sum(w)
// This is a product of leapfrog terms
// U(w) = e^{0.5*w*A} e^{w*B} r^{0.5*w*A}
// and for $n$ terms the expansion is U(w_1) U(w_2) ... U(w_{n-1}) U(w_n) U(w_{n-1}) ... U(w_1) ((m-1)/2 odd)
// or U(w_1) U(w_2) ... U(w_{n-1}) U(w_n) U(w_n) U(w_{n-1}) ... U(w_1) ((m-1)/2 even)
// We only need to supply (n-1) weights, as the final weight is 1 - 2*sum_{i<n}(w_i)
LTSDecomposition
LeapfrogDecompositionOdd(int Order, std::string Description, std::vector<double> w)
{
   std::vector<double> a;
   std::vector<double> b;
   //w.push_back(1.0 - 2.0 * std::accumulate(w.begin(), w.end(), 0.0));
   double aa = 0;
   for (auto ww : w)
   {
      a.push_back(0.5*(aa+ww));
      b.push_back(ww);
      aa = ww;
   }
   return SymmetricDecomposition(Order, Description, a, b);
}

std::map<std::string, LTSDecomposition>
Decompositions = {
   {"firstorder",     LTSDecomposition(1, "First order decomposition", {1.0}, {1.0})},
   {"leapfrog2",      LeapfrogDecompositionOdd(2, "Traditional 2nd order 3-term leapfrog decomposition", {})},
   {"optimized2-5",   SymmetricDecomposition(2, "Optimized 2nd order 5-term decomposition",
                                             {0.211324865405187}, {})},
   {"optimized4-11",  SymmetricDecomposition(4, "Optimized 4th order 11-term decomposition",
                                             {0.095848502741203681182, -0.078111158921637922695},
                                             {0.42652466131587616168, -0.12039526945509726545})},
   {"symmetric6-19",  LeapfrogDecompositionOdd(6, "6th order leapfrog 19-term decomposition",
                                               {0.18793069262651671457, 0.5553,
                                                     0.12837035888423653774, -0.84315275357471264676})}};


#if 0
void DoTEBD(StateComponent& A, StateComponent& B, RealDiagonalOperator& Lambda,
            SimpleOperator const& U, StatesInfo const& Info)
{
   // Algorithm avoids matrix inversion.  Let A,B be left-orthogonalized.
   // Let C^{s1,s2} = A^{s1} B^{s2}
   // Let C'^{s1,s2} = U(C^{s1,s2}) is the unitary gate applied to C

   // Let X^{s1,s2} = C'^{s1,s2} lambda_2

   // Do an SVD of X.  X^{s1,s2} = A'{s1} lambda_1 B'{s2}
   // Use orthogonality of A': sum_{s1} A'\dagger{s1} A'{s1} = I
   // sum_{s1} A'\dagger{s1} C'{s1,s2} = lambda_1 B'{s2} lambda_2^{-1} = G^{s2}

   // Then form D^{s2,s1} = G^{s2} A^{s1}
   // and so on

   StateComponent C = local_tensor_prod(A,B);
   StateComponent Cu = local_prod(U, C);

   StateComponent X = Cu * Lambda;
   AMatSVD SL(X, Tensor::ProductBasis<BasisList, BasisList>(A.LocalBasis(), B.LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                               SInfo, Info);
   std::cout << " Entropy=" << Info.TotalEntropy()
             << " States=" << Info.KeptStates()
             << " Trunc=" << Info.TruncationError()
             << std::endl;

   SL.ConstructMatrices(SL.begin(), Cutoff, A, Lambda, B);

   StateComponent G = partial_prod(herm(A), Cu,
                                   Tensor::ProductBasis<BasisList, BasisList>(A.LocalBasis(),
                                                                              B.LocalBasis()));

   B = A;
   A = Lambda*G;
}
#else
LinearAlgebra::DiagonalMatrix<double>
InvertDiagonal(LinearAlgebra::DiagonalMatrix<double> const& D, double Tol = 1E-15)
{
   LinearAlgebra::DiagonalMatrix<double> Result(size1(D), size2(D));
   for (unsigned i = 0; i < size1(D); ++i)
   {
      Result.diagonal()[i] = norm_frob(D.diagonal()[i]) < Tol ? 0.0 : 1.0 / D.diagonal()[i];
   }
   return Result;
}

RealDiagonalOperator
InvertDiagonal(RealDiagonalOperator const& D, double Tol = 1E-15)
{
   RealDiagonalOperator Result(D.Basis1(), D.Basis2());
   for (unsigned i = 0; i < D.Basis1().size(); ++i)
   {
      Result(i,i) = InvertDiagonal(D(i,i), Tol);
   }
   return Result;
}

//
// Input is A, B, Lambda
// A,B are in left-canonical form
// Lambda01 Gamma1 Lambda12 Gamma2 Lambda23
// A = Lambda01 Gamma1
// B = Lambda12 Gamma2
// Lambda = Lambda23
//
// On exit, A and B are in left canonical form,
// final Lambda' = Lambda12
//

TruncationInfo
DoTEBD(StateComponent& A, StateComponent& B, RealDiagonalOperator& Lambda,
       SimpleOperator const& U, StatesInfo const& SInfo)
{
   // simple algorithm with matrix inversion
   RealDiagonalOperator LambdaSave = Lambda;
   StateComponent C = local_tensor_prod(A,B);
   StateComponent Cu = local_prod(U, C);
   StateComponent X = Cu * Lambda;
   AMatSVD SL(X, Tensor::ProductBasis<BasisList, BasisList>(A.LocalBasis(), B.LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                               SInfo, Info);
   SL.ConstructMatrices(SL.begin(), Cutoff, A, Lambda, B);
   // normalize
   Lambda *= 1.0 / norm_frob(Lambda);

   StateComponent G = B * InvertDiagonal(LambdaSave, 1E-8);

   B = Lambda*G;

   CHECK_EQUAL(A.Basis2(), B.Basis1());
   CHECK_EQUAL(A.Basis2(), Lambda.Basis1());

   return Info;
}
#endif

