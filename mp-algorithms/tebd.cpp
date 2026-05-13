// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/tebd.cpp
//
// Copyright (C) 2020-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp-algorithms/tebd.h"
#include "mp-algorithms/time-dependent-mpo.h"
#include "common/statistics.h"
#include "mps/density.h"
#include <stdexcept>

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

// A selection of decompositions.
// firstorder is teh standard first order decomposition
// leapfrog2 is the standard second order decomposition, error prefactor is (3/2)^2 / 8 = 0.28125
// optimized2-5 is a 5-step second order decomposition from
// Barthel and Zhang [https://doi.org/10.1016/j.aop.2020.168165  https://arxiv.org/abs/1901.04974]
// error prefactor 0.069778
// optimized4-11 from Barthel ansd Zhang Eqn (30a), error prefactor 0.10509.  Note that there is
// another 4-order decomposition with a smaller error prefactor Eqn (32) from Barthel and Zhang,
// but it isn't much improvement and at the cost of 13 terms versus 11.
// symmetric6-19 is Eqn (35) from Barthel and Zhang, with error prefactor 0.17255, negligably higher
// than their best 6th-order decomposition but 19 terms versus 23.
//
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

// FIXME: The finite and periodic gate assembly below intentionally duplicates
// some traversal logic while finite/infinite MPO and lattice boundary semantics
// are still represented separately.  Consolidate this once those abstractions
// make edge versus wraparound behavior explicit.

SimpleOperator
TEBDBondIdentity(OperatorComponent const& Left, OperatorComponent const& Right)
{
   return tensor_prod(SimpleOperator::make_identity(Left.LocalBasis1()),
                      SimpleOperator::make_identity(Right.LocalBasis1()));
}

SimpleOperator
ExponentiateTEBDBond(std::complex<double> Factor,
                     SimpleOperator const& BondTerm,
                     OperatorComponent const& Left,
                     OperatorComponent const& Right)
{
   if (BondTerm.is_null())
      return TEBDBondIdentity(Left, Right);
   return Exponentiate(Factor * BondTerm);
}

BasicTriangularMPO
PrepareFiniteTEBDHamiltonian(BasicTriangularMPO HamMPO, int PsiSize)
{
   CHECK(!HamMPO.empty())("Finite TEBD Hamiltonian MPO is empty");
   CHECK(HamMPO.size() <= PsiSize)("Finite TEBD Hamiltonian is larger than the wavefunction")
      (HamMPO.size())(PsiSize);
   CHECK(PsiSize % HamMPO.size() == 0)
      ("Finite TEBD wavefunction size must be a multiple of the Hamiltonian unit cell")
      (PsiSize)(HamMPO.size());

   if (HamMPO.size() < PsiSize)
      HamMPO = repeat(HamMPO, PsiSize / HamMPO.size());
   return HamMPO;
}

std::vector<SimpleOperator>
FiniteTEBDBondHamiltonian(BasicTriangularMPO const& HamMPO)
{
   if (HamMPO.size() < 2)
      return {};

   std::vector<SimpleOperator> BondH(HamMPO.size()-1);
   auto SplitTerms = SplitMPO(HamMPO);
   for (int i = 0; i < SplitTerms.size(); ++i)
   {
      for (auto const& op : SplitTerms[i])
      {
         if (op.size() == 1)
         {
            // Split single-site operators over the two bonds, unless they are at the edge.
            if (i == 0)
            {
               SimpleRedOperator p = tensor_prod(op[0](0,0), HamMPO[i+1](0,0));
               BondH[i] += p.scalar();
            }
            else if (i == HamMPO.size()-1)
            {
               SimpleRedOperator p = tensor_prod(HamMPO[i-1](0,0), op[0](0,0));
               BondH[i-1] += p.scalar();
            }
            else
            {
               SimpleRedOperator p = tensor_prod(op[0](0,0), HamMPO[i+1](0,0));
               BondH[i] += 0.5 * p.scalar();
               SimpleRedOperator q = tensor_prod(HamMPO[i-1](0,0), op[0](0,0));
               BondH[i-1] += 0.5 * q.scalar();
            }
         }
         else if (op.size() == 2)
         {
            // Ignore terms that cross beyond the finite-system boundary.
            if (i < HamMPO.size()-1)
               BondH[i] += coarse_grain(op).scalar();
         }
         else
         {
            throw std::runtime_error("fatal: operator has support over more than 2 sites.");
         }
      }
   }

   return BondH;
}

std::vector<SimpleOperator>
AssembleFiniteTEBDEvenSlice(BasicTriangularMPO const& HamMPO,
                            std::complex<double> SliceTimestep)
{
   std::vector<SimpleOperator> BondH = FiniteTEBDBondHamiltonian(HamMPO);
   std::vector<SimpleOperator> Terms;
   for (int i = 0; i < BondH.size(); i += 2)
   {
      Terms.push_back(ExponentiateTEBDBond(-SliceTimestep*std::complex<double>(0.0, 1.0), BondH[i],
                                           HamMPO[i], HamMPO[i+1]));
   }
   return Terms;
}

std::vector<SimpleOperator>
AssembleFiniteTEBDOddSlice(BasicTriangularMPO const& HamMPO,
                           std::complex<double> SliceTimestep)
{
   std::vector<SimpleOperator> BondH = FiniteTEBDBondHamiltonian(HamMPO);
   std::vector<SimpleOperator> Terms;
   for (int i = 1; i < BondH.size(); i += 2)
   {
      Terms.push_back(ExponentiateTEBDBond(-SliceTimestep*std::complex<double>(0.0, 1.0), BondH[i],
                                           HamMPO[i], HamMPO[i+1]));
   }
   return Terms;
}

TEBDHamiltonianGates
AssembleFiniteTEBDHamiltonian(BasicTriangularMPO const& HamMPO,
                              std::complex<double> Timestep,
                              LTSDecomposition const& decomp)
{
   std::vector<std::vector<SimpleOperator>> EvenU;
   for (auto x : decomp.a())
      EvenU.push_back(AssembleFiniteTEBDEvenSlice(HamMPO, x*Timestep));

   std::vector<std::vector<SimpleOperator>> OddU;
   for (auto x : decomp.b())
      OddU.push_back(AssembleFiniteTEBDOddSlice(HamMPO, x*Timestep));

   std::vector<SimpleOperator> EvenContinuation;
   if (decomp.a().size() == decomp.b().size()+1)
      EvenContinuation = AssembleFiniteTEBDEvenSlice(HamMPO, (decomp.a().front() + decomp.a().back()) * Timestep);

   return {EvenU, OddU, EvenContinuation};
}

TEBDHamiltonianGates
AssembleFiniteTimeDependentTEBDHamiltonian(InfiniteLattice const& Lattice,
                                           std::string const& HamOperator,
                                           std::string const& TimeVar,
                                           std::complex<double> StepStart,
                                           std::complex<double> Timestep,
                                           LTSDecomposition const& decomp,
                                           int MagnusOrder,
                                           int MagnusQuadrature,
                                           int PsiSize)
{
   auto HamiltonianSlice = [&](std::complex<double> SliceStart, std::complex<double> SliceTimestep)
   {
      BasicTriangularMPO HamMPO = TimeDependentHamiltonianMPO(Lattice, HamOperator, TimeVar,
                                                              SliceStart, SliceTimestep,
                                                              MagnusOrder, MagnusQuadrature);
      return PrepareFiniteTEBDHamiltonian(HamMPO, PsiSize);
   };

   std::vector<std::vector<SimpleOperator>> EvenU;
   std::complex<double> EvenTime = StepStart;
   for (double x : decomp.a())
   {
      std::complex<double> SliceTimestep = x * Timestep;
      EvenU.push_back(AssembleFiniteTEBDEvenSlice(HamiltonianSlice(EvenTime, SliceTimestep), SliceTimestep));
      EvenTime += SliceTimestep;
   }

   std::vector<std::vector<SimpleOperator>> OddU;
   std::complex<double> OddTime = StepStart;
   for (double x : decomp.b())
   {
      std::complex<double> SliceTimestep = x * Timestep;
      OddU.push_back(AssembleFiniteTEBDOddSlice(HamiltonianSlice(OddTime, SliceTimestep), SliceTimestep));
      OddTime += SliceTimestep;
   }

   std::vector<SimpleOperator> EvenContinuation;
   if (decomp.a().size() == decomp.b().size()+1)
   {
      std::complex<double> SliceStart = StepStart - decomp.a().back() * Timestep;
      std::complex<double> SliceTimestep = (decomp.a().back() + decomp.a().front()) * Timestep;
      EvenContinuation = AssembleFiniteTEBDEvenSlice(HamiltonianSlice(SliceStart, SliceTimestep), SliceTimestep);
   }

   return {EvenU, OddU, EvenContinuation};
}

BasicTriangularMPO
PreparePeriodicTEBDHamiltonianUnitCell(BasicTriangularMPO HamMPO, int Coarsegrain, int PsiSize)
{
   HamMPO = coarse_grain(HamMPO, Coarsegrain);

   int UnitCellSize = statistics::lcm(PsiSize, HamMPO.size(), 2);
   return repeat(HamMPO, UnitCellSize / HamMPO.size());
}

std::vector<SimpleOperator>
PeriodicTEBDBondHamiltonian(BasicTriangularMPO const& HamMPO)
{
   int UnitCellSize = HamMPO.size();

   std::vector<SimpleOperator> BondH(UnitCellSize);
   auto SplitTerms = SplitMPO(HamMPO);
   for (int i = 0; i < SplitTerms.size(); ++i)
   {
      for (auto const& op : SplitTerms[i])
      {
         if (op.size() == 1)
         {
            // Split single-site operators over the two bonds, wrapping around the unit cell.
            SimpleRedOperator p = tensor_prod(op[0](0,0), HamMPO[(i+1)%UnitCellSize](0,0));
            BondH[i] += 0.5 * p.scalar();
            SimpleRedOperator q = tensor_prod(HamMPO[(i+UnitCellSize-1)%UnitCellSize](0,0), op[0](0,0));
            BondH[(i+UnitCellSize-1)%UnitCellSize] += 0.5 * q.scalar();
         }
         else if (op.size() == 2)
         {
            BondH[i] += coarse_grain(op).scalar();
         }
         else
         {
            throw std::runtime_error("fatal: operator has support over more than 2 sites.");
         }
      }
   }

   return BondH;
}

std::vector<SimpleOperator>
AssemblePeriodicTEBDEvenSlice(BasicTriangularMPO const& HamMPO,
                              std::complex<double> SliceTimestep)
{
   int UnitCellSize = HamMPO.size();
   std::vector<SimpleOperator> BondH = PeriodicTEBDBondHamiltonian(HamMPO);
   std::vector<SimpleOperator> Terms;
   for (int i = 0; i < BondH.size(); i += 2)
   {
      Terms.push_back(ExponentiateTEBDBond(-SliceTimestep*std::complex<double>(0.0, 1.0), BondH[i],
                                           HamMPO[i], HamMPO[(i+1)%UnitCellSize]));
   }
   return Terms;
}

std::vector<SimpleOperator>
AssemblePeriodicTEBDOddSlice(BasicTriangularMPO const& HamMPO,
                             std::complex<double> SliceTimestep)
{
   int UnitCellSize = HamMPO.size();
   std::vector<SimpleOperator> BondH = PeriodicTEBDBondHamiltonian(HamMPO);
   std::vector<SimpleOperator> Terms;

   // The odd slice uses a rotation of the unit cell that takes the right-most site and puts it on the left.
   Terms.push_back(ExponentiateTEBDBond(-SliceTimestep*std::complex<double>(0.0, 1.0), BondH[BondH.size()-1],
                                        HamMPO[BondH.size()-1], HamMPO[0]));
   for (int i = 1; i < BondH.size()-1; i += 2)
   {
      Terms.push_back(ExponentiateTEBDBond(-SliceTimestep*std::complex<double>(0.0, 1.0), BondH[i],
                                           HamMPO[i], HamMPO[(i+1)%UnitCellSize]));
   }

   return Terms;
}

TEBDHamiltonianGates
AssemblePeriodicTEBDHamiltonian(BasicTriangularMPO const& HamMPO,
                                std::complex<double> Timestep,
                                LTSDecomposition const& decomp)
{
   int UnitCellSize = HamMPO.size();
   std::vector<SimpleOperator> BondH = PeriodicTEBDBondHamiltonian(HamMPO);

   std::vector<std::vector<SimpleOperator>> EvenU;
   for (auto x : decomp.a())
   {
      std::vector<SimpleOperator> Terms;
      for (int i = 0; i < BondH.size(); i += 2)
      {
         Terms.push_back(ExponentiateTEBDBond(-Timestep*std::complex<double>(0,x), BondH[i],
                                              HamMPO[i], HamMPO[(i+1)%UnitCellSize]));
      }
      EvenU.push_back(std::move(Terms));
   }

   std::vector<std::vector<SimpleOperator>> OddU;
   for (auto x : decomp.b())
   {
      std::vector<SimpleOperator> Terms;
      Terms.push_back(ExponentiateTEBDBond(-Timestep*std::complex<double>(0,x), BondH[BondH.size()-1],
                                           HamMPO[BondH.size()-1], HamMPO[0]));
      for (int i = 1; i < BondH.size()-1; i += 2)
      {
         Terms.push_back(ExponentiateTEBDBond(-Timestep*std::complex<double>(0,x), BondH[i],
                                              HamMPO[i], HamMPO[(i+1)%UnitCellSize]));
      }
      OddU.push_back(std::move(Terms));
   }

   std::vector<SimpleOperator> EvenContinuation;
   if (decomp.a().size() == decomp.b().size()+1)
   {
      double x = decomp.a().front() + decomp.a().back();
      for (int i = 0; i < BondH.size(); i += 2)
      {
         EvenContinuation.push_back(ExponentiateTEBDBond(-Timestep*std::complex<double>(0,x), BondH[i],
                                                         HamMPO[i], HamMPO[(i+1)%UnitCellSize]));
      }
   }

   return {EvenU, OddU, EvenContinuation};
}

TEBDHamiltonianGates
AssemblePeriodicTimeDependentTEBDHamiltonian(InfiniteLattice const& Lattice,
                                             std::string const& HamOperator,
                                             std::string const& TimeVar,
                                             std::complex<double> StepStart,
                                             std::complex<double> Timestep,
                                             LTSDecomposition const& decomp,
                                             int MagnusOrder,
                                             int MagnusQuadrature,
                                             int Coarsegrain,
                                             int PsiSize)
{
   auto HamiltonianSlice = [&](std::complex<double> SliceStart, std::complex<double> SliceTimestep)
   {
      BasicTriangularMPO HamMPO = TimeDependentHamiltonianMPO(Lattice, HamOperator, TimeVar,
                                                              SliceStart, SliceTimestep,
                                                              MagnusOrder, MagnusQuadrature);
      return PreparePeriodicTEBDHamiltonianUnitCell(HamMPO, Coarsegrain, PsiSize);
   };

   std::vector<std::vector<SimpleOperator>> EvenU;
   std::complex<double> EvenTime = StepStart;
   for (double x : decomp.a())
   {
      std::complex<double> SliceTimestep = x * Timestep;
      EvenU.push_back(AssemblePeriodicTEBDEvenSlice(HamiltonianSlice(EvenTime, SliceTimestep), SliceTimestep));
      EvenTime += SliceTimestep;
   }

   std::vector<std::vector<SimpleOperator>> OddU;
   std::complex<double> OddTime = StepStart;
   for (double x : decomp.b())
   {
      std::complex<double> SliceTimestep = x * Timestep;
      OddU.push_back(AssemblePeriodicTEBDOddSlice(HamiltonianSlice(OddTime, SliceTimestep), SliceTimestep));
      OddTime += SliceTimestep;
   }

   std::vector<SimpleOperator> EvenContinuation;
   if (decomp.a().size() == decomp.b().size()+1)
   {
      std::complex<double> SliceStart = StepStart - decomp.a().back() * Timestep;
      std::complex<double> SliceTimestep = (decomp.a().back() + decomp.a().front()) * Timestep;
      EvenContinuation = AssemblePeriodicTEBDEvenSlice(HamiltonianSlice(SliceStart, SliceTimestep), SliceTimestep);
   }

   return {EvenU, OddU, EvenContinuation};
}

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
       double& LogAmplitude,
       SimpleOperator const& U, StatesInfo const& SInfo)
{
   // simple algorithm with matrix inversion
   Tensor::ProductBasis<BasisList, BasisList> PB(A.LocalBasis(), B.LocalBasis());
   RealDiagonalOperator LambdaSave = Lambda;
   StateComponent C = local_tensor_prod(A,B, PB);
   StateComponent Cu = local_prod(U, C);
   StateComponent X = Cu * Lambda;
   AMatSVD SL(X, PB);
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                               SInfo, Info);
   SL.ConstructMatrices(SL.begin(), Cutoff, A, Lambda, B);
   // normalize
   double a = norm_frob(Lambda);
   LogAmplitude += std::log(a);
   Lambda *= 1.0 / a;

   // todo: we can avoid the diagonal inversion here by constructing
   // +---A----
   // |   |
   // |   |
   // |   | |
   // +---Cu---
   // which is the left-orthogonal form of B.  This works because if you
   // imagine doing an SVD of Cu directly (without multiplying on the right by Lambda),
   // then both U and V are, in principle, left-orthogonalized.  The contraction over the
   // first part is the basis transformataion from this left A matrix to the basis
   // in which we are going to do the contraction.  This is the 'Hastings trick'.
   // StateComponent B = contract_local_tensor_prod_left(herm(A), Cu, PB);
   StateComponent G = B * InvertDiagonal(LambdaSave, 1E-8);
   B = Lambda*G;

   CHECK_EQUAL(A.Basis2(), B.Basis1());
   CHECK_EQUAL(A.Basis2(), Lambda.Basis1());

   return Info;
}
