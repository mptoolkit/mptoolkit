// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-iwigner-eckart.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

//
// Project a wavefunction using the Wigner-Eckart theorem
//
// The Regularize option is bugged - somehow the C_right and C_old
// matrices end up with incompatible bases, perhaps due to sorting of quantum numbers?

#include "wavefunction/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "tensor/wigner_eckart.h"
#include "tensor/regularize.h"

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;

StateComponent wigner_eckart(StateComponent const& A, SymmetryList const& FinalSL)
{
   // Get the projected local basis
   BasisList AbelianLocalBasis(FinalSL);
   for (unsigned i = 0; i < A.LocalBasis().size(); ++i)
   {
      QuantumNumbers::ProjectionList pl = enumerate_projections(A.LocalBasis()[i]);
      for (unsigned pi = 0; pi < pl.size(); ++pi)
      {
         AbelianLocalBasis.push_back(map_projection_to_quantum(pl[pi], FinalSL));
      }
   }

   StateComponent Result(AbelianLocalBasis, W1.AbelianBasis(), W2.AbelianBasis());
   int k = 0;
   for (unsigned i = 0; i < I->LocalBasis().size(); ++i)
   {
      QuantumNumbers::ProjectionList pl = enumerate_projections(I->LocalBasis()[i]);
      for (unsigned pi = 0; pi < pl.size(); ++pi)
      {
         Result[k++] = wigner_eckart((*I)[i], pl[pi], W1, W2);
      }
   }
   return Result;
}

InfiniteWavefunctionLeft
wigner_eckart(InfiniteWavefunctionLeft const& Psi, SymmetryList const& FinalSL)
{
   InfiniteWavefunctionLeft Result;

   VectorBasis b1 = Psi.Basis1();
   WignerEckartBasis<VectorBasis> W2(b1, FinalSL);

   Result.setBasis1(W2);

   // get the identity projection
   QuantumNumbers::ProjectionList PL = enumerate_projections(Psi.lambda_l().TransformsAs());
   CHECK_EQUAL(PL.size(), 1U);
   Projection IdentP = PL[0];

   QuantumNumbers::ProjectionList QPL = enumerate_projections(Psi.shift());
   Result.QShift = change(QuantumNumber(FinalSL), QPL[0]);

   WignerEckartBasis<VectorBasis> W2;
   for (int i = 0; i < Psi.size(); ++i)
   {

      StateComponent C = Psi[i];
      WignerEckartBasis<VectorBasis> W1 = W2;
      W2 = WignerEckartBasis<VectorBasis>(C.Basis2(), FinalSL);

      Result.push_back_lambda(wigner_eckart(Psi.lambda(i), IdentP, W1, W1));
      Result.push_back(wigner_eckart(C));
   }

   Result.push_back_lambda(wigner_eckart(Psi.lambda_r(), IdentP, W2, W2));

   Result.setBasis2(W2);

   Result.check_structure();

   return Result;
}

// functor to use the visitor pattern with wavefunction types
struct ApplyWignerEckart : public boost::static_visitor<WavefunctionTypes>
{
   ApplyWignerEckart(SymmetryList const& FinalSL_)
      : FinalSL(FinalSL_) {}

   template <typename T>
   T operator()(T const& Psi) const
   {
      return wigner_eckart(T, FinalSL);
   }

   SymmetryList FinalSL;
};

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-wigner-eckart <input-psi> <output-psi> <symmetry-list>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[2], 1, PageSize, CacheSize);
   pvalue_ptr<MPWavefunction> PsiIn = pheap::ImportHeap(argv[1]);
   std::string FinalSLStr = argv[3];

   SymmetryList FinalSL = SymmetryList(FinalSLStr);

   MPWavefunction Result;
   Result.AppendHistory(EscapeCommandline(argv, arc));
   Result.Wavefunction() = boost::apply_visitor(ApplyWignerEckart(FinalSL), PsiIn->Wavefunction());

   pvalue_ptr<MPWavefunction> PsiNew = new MPWavefunction(Result);

   pheap::ShutdownPersistent(PsiNew);
}
