// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/random_wavefunc.cpp
//
// Copyright (C) 2004-2023 Ian McCulloch <ian@qusim.net>
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

#include "random_wavefunc.h"
#include "linearalgebra/matrix_utility.h"
#include "common/randutil.h"

bool
WavefunctionDesc::Flip(std::vector<BasisList> const& Basis, int Site, int NewState,
                       QuantumNumber const& NewQuantumNumber)
{
   Projection Delta = difference(NewQuantumNumber, Height[Site]);
   //   TRACE(Delta);
   State[Site] = NewState;
   for (int i = Site; i >= 0; --i)
   {
      if (!is_possible(Height[i], Delta)) return false;
      //      TRACE(i)(Height[i])(change(Height[i], Delta));

      QuantumNumbers::QuantumNumberList N = transform_targets(Height[i], Basis[i][State[i]]);
      if (N.size() < 1) return false;

      Height[i] = change(Height[i], Delta);

      if (!is_transform_target(Height[i+1], Basis[i][State[i]], Height[i]))
      {
         //DEBUG_TRACE(Height[i+1])(Basis[i][State[i]])(Height[i]);
         return false;
      }
   }
   return true;
}

WavefunctionDesc::WavefunctionDesc(std::vector<BasisList> const& L)
   : State(L.size()), Height(L.size()+1)
{
   QuantumNumber Q(L[0].GetSymmetryList());
   Height[L.size()] = Q;
   for (int i = L.size()-1; i >= 0; --i)
   {
      BasisList Basis = L[i];
      State[i] = randutil::rand_int(0, Basis.size()-1);

      QuantumNumbers::QuantumNumberList QList = transform_targets(Q, Basis[State[i]]);
      Q = QList[randutil::rand_int(0, QList.size()-1)];
      Height[i] = Q;
   }
}

WavefunctionDesc::WavefunctionDesc(std::vector<BasisList> const& L,
                                   QuantumNumber Q)
   : State(L.size()), Height(L.size()+1)
{
   Height[L.size()] = Q;
   for (int i = L.size()-1; i >= 0; --i)
   {
      BasisList Basis = L[i];
      State[i] = randutil::rand_int(0, Basis.size()-1);

      QuantumNumbers::QuantumNumberList QList = transform_targets(Q, Basis[State[i]]);
      Q = QList[randutil::rand_int(0, QList.size()-1)];
      Height[i] = Q;
   }
}

std::ostream& operator<<(std::ostream& out, WavefunctionDesc const& Config)
{
   int i = Config.State.size();
   out << "Height: " << Config.Height[i] << '\n';
   while (i > 0)
   {
      --i;
      out << "State: " << Config.State[i] << '\n';
      out << "Height: " << Config.Height[i] << '\n';
   }
   return out;
}

void WavefunctionDesc::CheckValid(std::vector<BasisList> const& Basis) const
{
   for (int i = Basis.size()-1; i >= 0; --i)
   {
      CHECK(is_projection(Basis[i][State[i]], difference(Height[i], Height[i+1])));
   }
}

WavefunctionDesc
CreateRandomConfiguration(std::vector<BasisList> const& Basis,
                          QuantumNumber const& q, double Beta,
                          QuantumNumber RightBoundary)
{
   if (RightBoundary.is_null())
      RightBoundary = QuantumNumber(Basis[0].GetSymmetryList());
   double now = clock() / CLOCKS_PER_SEC;
   WavefunctionDesc Psi(Basis, RightBoundary);
   while (q != Psi.TransformsAs())
   {
      double c = weight(difference(q, Psi.TransformsAs()));
      int Site = randutil::rand_int(0, Basis.size()-1);
      int NewState = randutil::rand_int(0, Basis[Site].size()-1);
      QuantumNumbers::QuantumNumberList QList = transform_targets(Psi.Height[Site+1], Basis[Site][NewState]);
      //      for (std::size_t Qi = 0; Qi < QList.size(); ++Qi)
      {
         QuantumNumber NewQ(QList[randutil::rand_int(0, QList.size()-1)]);
         WavefunctionDesc New = Psi;
         if (New.Flip(Basis, Site, NewState, NewQ))
         {
            //      TRACE(Psi.TransformsAs())(New.TransformsAs());
            if (clock() / CLOCKS_PER_SEC > now+10)
            {
               now = clock() / CLOCKS_PER_SEC;
               TRACE("slow convergence!")(Site)(q)(NewQ)(New.TransformsAs())(Psi.TransformsAs());
            }
            double Newc = weight(difference(q, New.TransformsAs()));
            //TRACE(Newc);
            if (randutil::rand() < exp(Beta * (c-Newc)))
            {
               //TRACE("Accepted");
               Psi = New;
               //              continue;
            }
         }
      }
   }
   return Psi;
}

LinearWavefunction WavefunctionFromConfiguration(WavefunctionDesc const& Psi, std::vector<BasisList> const& Basis,
						 QuantumNumber const& RightBoundary)
{
   BasisList Vac = make_single_basis(RightBoundary);
   VectorBasis B2(Vac);
   QuantumNumber Ident(RightBoundary.GetSymmetryList());  // the scalar quantum number

   LinearWavefunction Result;

   MatrixOperator Center(B2, B2, Ident);
   Center(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
   for (int i = Basis.size()-1; i >= 0; --i)
   {
      //      TRACE(i)(Psi.Height[i]);
      VectorBasis B1(Basis[i].GetSymmetryList());
      B1.push_back(Psi.Height[i], 1);
      MatrixOperator Next(B1, B2, Basis[i][Psi.State[i]]);
      Next(0,0) = LinearAlgebra::Matrix<double>(1,1,1);
      StateComponent R(Basis[i], B1, B2);
      R[Psi.State[i]] = Next;
      R = prod(R, Center);
      Center = Multiply(TruncateBasis1(R));
      Result.push_front(R);
      B2 = B1;
   }
   return Result;
}

LinearWavefunction CreateRandomWavefunction(std::vector<BasisList> const& Basis,
					    QuantumNumber const& q, double Beta,
					    QuantumNumber const& RightBoundary, int NConfig, int Verbose)
{
   LinearWavefunction Result;
   for (int i = 0; i < NConfig; ++i)
   {
      WavefunctionDesc Psi = CreateRandomConfiguration(Basis, q, Beta, RightBoundary);
      Result = Result + (randutil::rand()*2-1) * WavefunctionFromConfiguration(Psi, Basis, RightBoundary);
      if (Verbose > 0)
	 std::cout << "." << std::flush;
   }
   return Result;
}

LinearWavefunction CreateRandomWavefunction(std::vector<BasisList> const& Basis,
                                        QuantumNumber const& q, double Beta)
{
   return CreateRandomWavefunction(Basis, q, Beta, QuantumNumber(Basis[0].GetSymmetryList()));
}
