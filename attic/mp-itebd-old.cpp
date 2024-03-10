// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/mp-itebd-old.cpp
//
// Copyright (C) 2004-2021 Ian McCulloch <ian@qusim.net>
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

#include "mpo/basic_triangular_mpo.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/prog_options.h"
#include "lattice/unitcell_mpo.h"
#include "lattice/unitcell-parser.h"
#include <cctype>
#include "interface/inittemp.h"

namespace prog_opt = boost::program_options;

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

// generic Lie-Trotter-Suzuki decomposition with two slices, A and B.
class LTSDecomposition
{
   public:
      LTSDecomposition() : Order_(0) {}
      LTSDecomposition(int Order, std::string Description,
                       std::vector<double> a,
                       std::vector<double> b)
         : Order_(Order), Description_(Description), a_(a), b_(b)
      {
         CHECK(a.size() == b.size() || a.size() == b.size()+1);
         CHECK_CLOSE(std::accumulate(a.begin(), a.end(), 0.0), 1.0);
         CHECK_CLOSE(std::accumulate(b.begin(), b.end(), 0.0), 1.0);
      }

      int Order_;
      std::string Description_;
      std::vector<double> a_;
      std::vector<double> b_;
};

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

void DoEvenSlice(std::deque<StateComponent>& Psi,
                 std::deque<RealDiagonalOperator>& Lambda,
                 SimpleOperator const& UEven,
                 StatesInfo const& SInfo,
                 int Verbose)
{
   // physical state is A B (Lambda) C D (Lambda) ...
   // All A,B,C,D,... are left canonical
   // After even slice, the state is
   // A (Lambda) B C (Lambda) D E (Lambda) ...
   unsigned Sz = Psi.size();
   int MaxStates = 0;
   for (unsigned i = 0; i < Sz; i += 2)
   {
      TruncationInfo Info = DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UEven, SInfo);
      if (Verbose > 0)
      {
         std::cout << "Bond=" << (i+1)
                   << " Entropy=" << Info.TotalEntropy()
                   << " States=" << Info.KeptStates()
                   << " Trunc=" << Info.TruncationError()
                   << '\n';
      }
      MaxStates = std::max(MaxStates, Info.KeptStates());
   }
   std::cout << "Even slice finished, max states=" << MaxStates << '\n';
   std::cout << std::flush;
}

void DoOddSlice(std::deque<StateComponent>& Psi,
                std::deque<RealDiagonalOperator>& Lambda,
                QuantumNumber const& QShift,
                SimpleOperator const& UOdd,
                StatesInfo const& SInfo,
                int Verbose)
{
   // In preparation for the odd slice, we rotate to
   // Z A (Lambda) B C (Lambda) D E (Lambda) ...
   // After the odd slice, we have
   // Z (Lamda) A B (Lambda) C D (Lambda) ...
   // which we need to rotate back to
   // A B (Lambda) C D (Lambda) ..... Z (Lambda)
   unsigned Sz = Psi.size();
   Psi.push_front(delta_shift(Psi.back(), QShift));
   Psi.pop_back();
  int MaxStates = 0;
   for (unsigned i = 0; i < Sz; i += 2)
   {
      TruncationInfo Info = DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UOdd, SInfo);
      if (Verbose > 0)
      {
         std::cout << "Bond=" << i
                   << " Entropy=" << Info.TotalEntropy()
                   << " States=" << Info.KeptStates()
                   << " Trunc=" << Info.TruncationError()
                   << '\n';
      }
      MaxStates = std::max(MaxStates, Info.KeptStates());
   }
   std::cout << "Odd slice finished, max states=" << MaxStates << '\n';
   std::cout << std::flush;

   Psi.push_back(delta_shift(Psi.front(), adjoint(QShift)));
   Psi.pop_front();
   Lambda.push_back(delta_shift(Lambda.front(), adjoint(QShift)));
   Lambda.pop_front();
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string TimestepStr;  // so we can get the number of digits to use
      std::string InitialTimeStr;
      int N = 1;
      int SaveEvery = 1;
      std::string OpStr;
      std::string InputFile;
      std::string OutputPrefix;
      std::string Operator;
      std::string Operator2;
      int MinStates = 1;
      int States = 100000;
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      int OutputDigits = 0;
      int Coarsegrain = 1;
      std::string DecompositionStr = "leapfrog2";

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("wavefunction,w", prog_opt::value(&InputFile), "input wavefunction")
	 ("output,o", prog_opt::value(&OutputPrefix), "prefix for saving output files")
	 ("operator", prog_opt::value(&Operator), "operator for the first slice ST decomposition")
	 ("operator2", prog_opt::value(&Operator2), "operator for the second slice ST decomposition [optional]")
	 ("timestep,t", prog_opt::value(&TimestepStr), "timestep (required)")
         ("decomposition,c", prog_opt::value(&DecompositionStr), FormatDefault("choice of decomposition", DecompositionStr).c_str())
	 ("num-timesteps,n", prog_opt::value(&N), FormatDefault("number of timesteps to calculate", N).c_str())
	 ("save-timesteps,s", prog_opt::value(&SaveEvery), "save the wavefunction every s timesteps")
         ("min-states", prog_opt::value<int>(&MinStates),
          FormatDefault("Minimum number of states to keep", MinStates).c_str())
         ("states,m", prog_opt::value(&States),
          FormatDefault("number of states", States).c_str())
         ("trunc,r", prog_opt::value<double>(&TruncCutoff),
          FormatDefault("Truncation error cutoff", TruncCutoff).c_str())
         ("eigen-cutoff,d", prog_opt::value(&EigenCutoff),
          FormatDefault("Cutoff threshold for density matrix eigenvalues", EigenCutoff).c_str())
         ("coarsegrain", prog_opt::value(&Coarsegrain),
          "coarse-grain N-to-1 sites")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1
          || vm.count("operator") < 1 || vm.count("timestep") < 1)
      {
         print_copyright(std::cerr, "tools", "mp-itebd");
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "\nAvailable decompositions:\n";
         for (auto const& d : Decompositions)
         {
            std::cerr << d.first << " : " << d.second.Description_ << '\n';
         }
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      LTSDecomposition decomp;
      for (auto const& d : Decompositions)
      {
         if (d.first == DecompositionStr)
            decomp = d.second;
      }
      if (decomp.Order_ == 0)
      {
         std::cerr << "mp-itebd: fatal: invalid decomposition\n";
         return 1;
      }

      if (Verbose > 0)
         std::cout << "Loading wavefunction..." << std::endl;

      mp_pheap::InitializeTempPHeap();
      pvalue_ptr<MPWavefunction> PsiPtr = pheap::ImportHeap(InputFile);

      InfiniteWavefunctionLeft Psi = PsiPtr->get<InfiniteWavefunctionLeft>();

      if (OutputPrefix.empty())
         OutputPrefix = PsiPtr->Attributes()["Prefix"].as<std::string>();

      if (OutputPrefix.empty())
      {
         std::cerr << "mp-itebd: fatal: no output prefix specified\n";
         return 1;
      }

      if (InitialTimeStr.empty())
         InitialTimeStr = PsiPtr->Attributes()["Time"].as<std::string>();

      double InitialTime = InitialTimeStr.empty() ? 0.0 : stod(InitialTimeStr);

      double Timestep = stod(TimestepStr);

      if (OutputDigits == 0)
      {
         OutputDigits = std::max(Digits(InitialTimeStr), Digits(TimestepStr));
      }

      UnitCellMPO EvenOp, OddOp;
      InfiniteLattice Lattice;
      std::tie(EvenOp, Lattice) = ParseUnitCellOperatorAndLattice(Operator);
      if (Operator2.empty())
      {
	 OddOp = translate(EvenOp, Coarsegrain);
      }
      else
      {
	 std::tie(OddOp, Lattice) = ParseUnitCellOperatorAndLattice(Operator2);
      }

      EvenOp.ExtendToCover(2*Coarsegrain, 0);
      OddOp.ExtendToCover(2*Coarsegrain, Coarsegrain);

      if (EvenOp.offset() != 0 || EvenOp.size() != 2*Coarsegrain)
      {
	 std::cerr << "mp-itebd: fatal: slice 1 operator not valid.\n";
	 return 1;
      }

      if (OddOp.offset() != Coarsegrain || OddOp.size() != 2*Coarsegrain)
      {
	 std::cerr << "mp-itebd: fatal: slice 2 operator not valid.\n";
	 return 1;
      }

      SimpleOperator EvenX = coarse_grain(EvenOp.MPO()).scalar();
      SimpleOperator OddX = coarse_grain(EvenOp.MPO()).scalar();

      std::vector<SimpleOperator> EvenU;
      for (auto x : decomp.a_)
      {
         EvenU.push_back(Exponentiate(std::complex<double>(0, -Timestep * x) * EvenX));
      }

      std::vector<SimpleOperator> OddU;
      for (auto x : decomp.b_)
      {
         OddU.push_back(Exponentiate(std::complex<double>(0, -Timestep * x) * OddX));
      }

      // If we have an odd number of terms (the usual case), then we can wrap around the last even slice
      // if we are continuing the evolution beyond the current timestep.  This is the sum of a[last] + a[0] terms
      SimpleOperator EvenContinuation;
      if (decomp.a_.size() == decomp.b_.size()+1)
      {
         double x = decomp.a_.front() + decomp.a_.back();
         EvenContinuation = Exponentiate(std::complex<double>(0, -Timestep * x) * EvenX);
      }

      std::cout << "Using decomposition " << DecompositionStr << '\n';
      if (Verbose > 0)
      {
         std::cout << "Number of even slices: " << EvenU.size() << '\n';
         std::cout << "Number of odd slices: " << OddU.size() << '\n';
      }

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = States;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      std::cout << SInfo << '\n';

      if (Psi.size()%2 != 0)
      {
         std::cerr << "mp-itebd: warning: wavefunction is not a multiple of 2 sites, doubling unit cell...\n";
         Psi = repeat(Psi, 2);
      }

      QuantumNumber QShift = Psi.qshift();

      std::deque<StateComponent> PsiVec(Psi.begin(), Psi.end());
      std::deque<RealDiagonalOperator> Lambda;

      for (int i = 2; i <= Psi.size(); i += 2)
      {
         Lambda.push_back(Psi.lambda(i));
      }

      if (SaveEvery == 0)
         SaveEvery = N;

      // If Continue then we merge the final (even) slice of one timestep with the first slice of
      // the next timestep
      bool Continue = false;
      int tstep = 0;

      while (tstep < N)
      {
         if (Continue)
         {
            if (Verbose > 1)
            {
               std::cout << "Merge slice\n";
            }
            DoEvenSlice(PsiVec, Lambda, EvenContinuation, SInfo, Verbose);
         }
         else
         {
            DoEvenSlice(PsiVec, Lambda, EvenU[0], SInfo, Verbose);
         }

         DoOddSlice(PsiVec, Lambda, QShift, OddU[0], SInfo, Verbose);
         for (int bi = 1; bi < OddU.size(); ++bi)
         {
            DoEvenSlice(PsiVec, Lambda, EvenU[bi], SInfo, Verbose);
            DoOddSlice(PsiVec, Lambda, QShift, OddU[bi], SInfo, Verbose);
         }

         // do we do a continuation?
         Continue = EvenU.size() > OddU.size();

         ++tstep;
         std::cout << "Timestep " << tstep << " time " << (InitialTime+tstep*Timestep) << '\n';

         // do we save the wavefunction?
         if ((tstep % SaveEvery) == 0 || tstep == N)
         {
            LinearWavefunction Psi;
            if (EvenU.size() > OddU.size())
            {
               CHECK_EQUAL(EvenU.size(), OddU.size()+1);
               std::cout << "Doing final slice before saving wavefunction.\n";
               // do the final slice to finish the timstep.  Make a copy of the wavefunction since it is better to
               // avoid a truncation step and 'continue' the wavefunction by wrapping around the next timestep.
               std::deque<StateComponent> PsiVecSave = PsiVec;
               std::deque<RealDiagonalOperator> LambdaSave = Lambda;
               DoEvenSlice(PsiVecSave, LambdaSave, EvenU.back(), SInfo, Verbose);
               Psi = LinearWavefunction::FromContainer(PsiVecSave.begin(), PsiVecSave.end());
            }
            else
               Psi = LinearWavefunction::FromContainer(PsiVec.begin(), PsiVec.end());

            // save the wavefunction
            std::cout << "Saving wavefunction\n";
            MPWavefunction Wavefunction;
            std::string TimeStr = FormatDigits(InitialTime + tstep * Timestep, OutputDigits);
            InfiniteWavefunctionLeft PsiL = InfiniteWavefunctionLeft::Construct(Psi, QShift);
            Wavefunction.Wavefunction() = std::move(PsiL);
            Wavefunction.AppendHistoryCommand(EscapeCommandline(argc, argv));
            Wavefunction.SetDefaultAttributes();
            Wavefunction.Attributes()["Time"] = TimeStr;
            Wavefunction.Attributes()["Prefix"] = OutputPrefix;
            *PsiPtr.mutate() = std::move(Wavefunction);
            pheap::ExportHeap(OutputPrefix + TimeStr, PsiPtr);
         }
      }
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      if (e.Why == "File exists")
         std::cerr << "Note: use --force (-f) option to overwrite.\n";
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      pheap::Cleanup();
      return 1;
   }
}
