// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//----------------------------------------------------------------------------

// Test-only fixture writer. Deliberately bypasses canonical-form invariants;
// it is not an installed toolkit command.
#include "wavefunction/mpwavefunction.h"
#include "wavefunction/infinitewavefunctionleft.h"
#include "wavefunction/infinitewavefunctionright.h"
#include "interface/inittemp.h"
#include <cmath>
#include <iostream>
#include <limits>

struct MutableFixture : InfiniteWavefunctionLeft
{
   explicit MutableFixture(InfiniteWavefunctionLeft const& Psi)
      : InfiniteWavefunctionLeft(Psi) {}

   using CanonicalWavefunctionBase::set;
   using CanonicalWavefunctionBase::set_lambda;
};

InfiniteWavefunctionLeft entangled_fixture()
{
   SymmetryList SL("Null:Null");
   QuantumNumber Q(SL);
   BasisList Local(SL);
   for (int i = 0; i < 4; ++i)
      Local.push_back(Q);
   VectorBasis Bond(SL);
   Bond.push_back(Q, 2);
   StateComponent A(Local, Bond, Bond);
   // A_ij = |i><j|/sqrt(2): the transfer map sends R to trace(R) I/2.
   for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
      {
         LinearAlgebra::Matrix<std::complex<double>> M(2, 2, 0.0);
         M(i,j) = 1.0 / std::sqrt(2.0);
         A[2*i+j](0,0) = M;
      }
   LinearWavefunction Psi(SL);
   Psi.push_back(A);
   Psi.push_back(A);
   RealDiagonalOperator Lambda = RealDiagonalOperator::make_identity(Bond);
   Lambda *= 1.0 / std::sqrt(2.0);
   return InfiniteWavefunctionLeft::ConstructFromOrthogonal(Psi, Q, Lambda);
}

InfiniteWavefunctionLeft symmetry_fixture(bool SU2)
{
   SymmetryList SL(SU2 ? "S:SU(2)" : "Sz:U(1)");
   QuantumNumber Zero(SL), Half(SL, "1/2");
   BasisList Local(SL);
   Local.push_back(SU2 ? Zero : Half);
   VectorBasis Left(SL), Right(SL);
   Left.push_back(Half, 1);
   Right.push_back(SU2 ? Half : Zero, 1);
   StateComponent A(Local, Left, Right);
   A[0](0,0) = LinearAlgebra::Matrix<std::complex<double>>(1, 1, 1.0);
   LinearWavefunction Psi(SL);
   Psi.push_back(A);
   RealDiagonalOperator Lambda = RealDiagonalOperator::make_identity(Right);
   Lambda *= 1.0 / std::sqrt(double(Right.total_degree()));
   return InfiniteWavefunctionLeft::ConstructFromOrthogonal(Psi, SU2 ? Zero : Half, Lambda);
}

template <typename T>
void save(std::string const& Name, T const& Psi)
{
   MPWavefunction Wavefunction;
   Wavefunction.Wavefunction() = Psi;
   Wavefunction.SetDefaultAttributes();
   pheap::ExportObject(Name, Wavefunction);
}

int main()
{
   try
   {
      mp_pheap::InitializeTempPHeap();
      auto Good = entangled_fixture();
      save("canonical", Good);

      MutableFixture Bad(Good);
      Bad.set(0, 2.0 * Good[0]);
      save("tensor-error", static_cast<InfiniteWavefunctionLeft const&>(Bad));

      Bad = MutableFixture(Good);
      auto InteriorLambda = Good.lambda(1);
      InteriorLambda *= 2.0;
      Bad.set_lambda(1, InteriorLambda);
      save("stored-density-error", static_cast<InfiniteWavefunctionLeft const&>(Bad));

      Bad = MutableFixture(Good);
      auto Lambda = Good.lambda_r();
      Lambda(0,0)(0,0) = std::sqrt(0.9);
      Lambda(0,0)(1,1) = std::sqrt(0.1);
      Bad.set_lambda(0, Lambda);
      Bad.set_lambda(Good.size(), Lambda);
      save("fixed-point-error", static_cast<InfiniteWavefunctionLeft const&>(Bad));

      Bad = MutableFixture(Good);
      InteriorLambda *= 0.0;
      Bad.set_lambda(1, InteriorLambda);
      save("zero-reference", static_cast<InfiniteWavefunctionLeft const&>(Bad));

      Bad = MutableFixture(Good);
      Bad.set(0, std::numeric_limits<double>::quiet_NaN() * Good[0]);
      save("nonfinite", static_cast<InfiniteWavefunctionLeft const&>(Bad));

      Good.scale(std::complex<double>(0.0, 3.0));
      save("scaled", Good);
      save("right-canonical", InfiniteWavefunctionRight(Good));
      save("u1-shift", symmetry_fixture(false));
      auto SU2 = symmetry_fixture(true);
      save("su2", SU2);
      Bad = MutableFixture(SU2);
      Bad.set(0, 2.0 * SU2[0]);
      save("su2-tensor-error", static_cast<InfiniteWavefunctionLeft const&>(Bad));
      pheap::Shutdown();
   }
   catch (std::exception const& e)
   {
      std::cerr << e.what() << '\n';
      pheap::Cleanup();
      return 1;
   }
}
