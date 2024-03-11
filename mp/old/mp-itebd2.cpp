// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/old/mp-itebd2.cpp
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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

void DoTEBD(std::deque<StateComponent>& L, std::deque<StateComponent>& R, RealDiagonalOperator& Lambda,
            QuantumNumber const& QShift, int Polarity,
            SimpleOperator const& U, StatesInfo const& SInfo)
{
   StateComponent A = L.back(); L.pop_back();
   StateComponent B = R.front(); R.pop_front();
   DoTEBD(A, B, Lambda, QShift, Polarity, U, SInfo);

   std::swap(L, R);
   R.push_back(A);
   L.push_front(B);
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

void DoTEBD(StateComponent& A, StateComponent& B, RealDiagonalOperator& Lambda,
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
   std::cout << " Entropy=" << Info.TotalEntropy()
             << " States=" << Info.KeptStates()
             << " Trunc=" << Info.TruncationError()
             << std::endl;

   SL.ConstructMatrices(SL.begin(), Cutoff, A, Lambda, B);

   StateComponent G = B * InvertDiagonal(LambdaSave, 1E-8);

   B = Lambda*G;

   CHECK_EQUAL(A.Basis2(), B.Basis1());
   CHECK_EQUAL(A.Basis2(), Lambda.Basis1());
}
#endif

LinearWavefunction
DoEvenStep(std::vector<StateComponent> Psi, std::deque<RealDiagonalOperator> Lambda,
           SimpleOperator const& UEven, StatesInfo const& SInfo)
{
   for (unsigned i = 0; i < Psi.size(); i += 2)
   {
      DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UEven, SInfo);
   }
   return LinearWavefunction(Psi.begin(), Psi.end());
}

void DoEvenOddStep(std::deque<StateComponent>& Psi, std::deque<RealDiagonalOperator>& Lambda,
                   QuantumNumber const& QShift,
                   SimpleOperator const& UEven, SimpleOperator const& UOdd, StatesInfo const& SInfo)
{
   // physical state is A B (Lambda) C D (Lambda) ...
   // All A,B,C,D,... are left canonical
   // After even slice, the state is
   // A (Lambda) B C (Lambda) D E (Lambda) ...
   // In preparation for the odd slice, we rotate to
   // Z A (Lambda) B C (Lambda) D E (Lambda) ...
   // After the odd slice, we have
   // Z (Lamda) A B (Lambda) C D (Lambda) ...
   // which we need to rotate back to
   // A B (Lambda) C D (Lambda) ..... Z (Lambda)
   for (unsigned i = 0; i < Psi.size(); i += 2)
   {
      DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UEven, SInfo);
   }
   // need to the pair that wraps around separately
   Psi.push_front(delta_shift(Psi.back(), QShift));
   Psi.pop_back();
   for (unsigned i = 0; i < Psi.size(); i += 2)
   {
      DoTEBD(Psi[i], Psi[i+1], Lambda[i/2], UOdd, SInfo);
   }
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
      int States = 100;
      double TruncCutoff = 0;
      double EigenCutoff = 1E-16;
      int OutputDigits = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("wavefunction,w", prog_opt::value(&InputFile), "input wavefunction")
	 ("output,o", prog_opt::value(&OutputPrefix), "prefix for saving output files")
	 ("operator", prog_opt::value(&Operator), "operator for the first slice ST decomposition")
	 ("operator2", prog_opt::value(&Operator2), "operator for the second slice ST decomposition [optional]")
	 ("timestep,t", prog_opt::value(&TimestepStr), "timestep (required)")
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
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

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
	 OddOp = translate(EvenOp, 1);
      }
      else
      {
	 std::tie(EvenOp, Lattice) = ParseUnitCellOperatorAndLattice(Operator2);
      }

      EvenOp.ExtendToCover(2, 0);
      OddOp.ExtendToCover(2, 1);

      if (EvenOp.offset() != 0 || EvenOp.size() != 2)
      {
	 std::cerr << "mp-itebd: fatal: slice 1 operator not valid.\n";
	 return 1;
      }

      if (OddOp.offset() != 1 || OddOp.size() != 2)
      {
	 std::cerr << "mp-itebd: fatal: slice 2 operator not valid.\n";
	 return 1;
      }

      SimpleOperator EvenX = coarse_grain(EvenOp.MPO()).scalar();
      SimpleOperator OddX = coarse_grain(EvenOp.MPO()).scalar();

      SimpleOperator EvenU = Exponentiate(std::complex<double>(0, -Timestep) * EvenX);
      SimpleOperator OddU = Exponentiate(std::complex<double>(0, -Timestep) * OddX);

      SimpleOperator EvenUHalf = Exponentiate(std::complex<double>(0, -Timestep/2) * EvenX);

      StatesInfo SInfo;
      SInfo.MinStates = MinStates;
      SInfo.MaxStates = States;
      SInfo.TruncationCutoff = TruncCutoff;
      SInfo.EigenvalueCutoff = EigenCutoff;

      std::cout << SInfo << '\n';

      if (Psi.size() == 1)
         Psi = repeat(Psi, 2);

      if (Psi.size()%2 != 0)
      {
         std::cerr << "mp-itebd: fatal: wavefunction must be multiple of 2 sites.\n";
         return 1;
      }

      QuantumNumber QShift = Psi.qshift();


      std::vector<StateComponent> PsiVec(Psi.begin(), Psi.end());
      std::deque<RealDiagonalOperator> Lambda;

      for (unsigned i = 2; i <= Psi.size(); ++i)
      {
         Lambda.push_back(Psi.lambda(2*(i+1)));
      }

      if (SaveEvery == 0)
         SaveEvery = N;

      // the initial half timestep
      int tstep = 1;
      DoEvenOddStep(Psi, Lambda, QShift, EvenUHalf, OddU, SInfo);
      std::cout << "Timestep " << tstep << " time " << (InitialTime+tstep*Timestep) << '\n';

      while (tstep < N)
      {
         while (tstep % SaveEvery != 0)
         {
            DoEvenOddStep(Psi, Lambda, QShift, EvenU, OddU, SInfo);
            ++tstep;
            std::cout << "Timestep " << tstep << " time " << (InitialTime+tstep*Timestep) << '\n';
         }

         // evolve the final half step and save the wavefunction
         // do this with a temporary copy, since if we continue then we do another slice anyway
         StateComponent Ax = A;
         StateComponent Bx = B;
         RealDiagonalOperator Lambdax = Lambda;

         DoTEBD(Ax, Bx, Lambdax, QShift, 1, EvenUHalf, SInfo);

         // save the wavefunction
         std::cout << "Saving wavefunction\n";
         LinearWavefunction Psi;
         Psi.push_back(Ax);
         Psi.push_back(Bx);
         MPWavefunction Wavefunction;
         std::string TimeStr = FormatDigits(InitialTime + tstep * Timestep, OutputDigits);
         InfiniteWavefunctionLeft PsiL = InfiniteWavefunctionLeft::Construct(Psi, QShift);
         // rotate to the A,B back into the right order
         PsiL.rotate_left(1);
         Wavefunction.Wavefunction() = std::move(PsiL);
         Wavefunction.AppendHistory(EscapeCommandline(argc, argv));
         Wavefunction.SetDefaultAttributes();
         Wavefunction.Attributes()["Time"] = TimeStr;
         Wavefunction.Attributes()["Prefix"] = OutputPrefix;
         *PsiPtr.mutate() = std::move(Wavefunction);

         pheap::ExportHeap(OutputPrefix + TimeStr, PsiPtr);

         if (tstep+1 < N)
         {
            DoTEBD(A, B, Lambda, QShift, 1, EvenU, SInfo);
            DoTEBD(A, B, Lambda, QShift, -1, OddU, SInfo);
            ++tstep;
            std::cout << "Timestep " << tstep << " time " << (InitialTime+tstep*Timestep) << '\n';
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
