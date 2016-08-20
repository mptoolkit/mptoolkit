// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/local-evolution.cpp
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

#include "local-evolution.h"

bool ShowStates = true;
bool ShowEntropy = true;
bool ShowTruncation = true;
bool ShowSmallestKept = true;

// assuming the input wavefunction is in normal form, applies a list of evolution operators
// to successive bonds.  For a suzuki-trotter decomposition, half of the evolution operators
// will the the identity operation - this should be negligible loss of efficiency as
// the actual application of the evolution term is O(d^2 m^2), versus the O(d^3 m^3)
// singular value decomposition.
void SweepRightEvolve(LinearWavefunction& Psi, std::list<SimpleOperator> const& BondOperators,
                      StatesInfo const& SInfo, bool ShowInfo)
{
   LinearWavefunction::iterator I = Psi.begin();
   LinearWavefunction::iterator J = I; ++J;
   std::list<SimpleOperator>::const_iterator BondIter = BondOperators.begin();
   MPStateComponent R = *I;

   int CurrentBond = 0;
   while (J != Psi.end())
   {
      CHECK(BondIter != BondOperators.end());
      MPStateComponent A = local_prod(*BondIter, local_tensor_prod(R, *J));
      AMatSVD SL(A, Tensor::ProductBasis<BasisList, BasisList>(R.LocalBasis(), J->LocalBasis()));
      TruncationInfo Info;
      AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                                  SInfo, Info);
      if (ShowInfo)
      {
         std::cout << "Bond=(" << (CurrentBond+1) << ',' << (CurrentBond+2) << ")";
         if (ShowStates) std::cout << " states=" << Info.KeptStates();
         if (ShowEntropy) std::cout << " entropy=" << Info.KeptEntropy();
         if (ShowTruncation) std::cout << " trunc=" << Info.TruncationError();
         if (ShowSmallestKept) std::cout << " skeep=" << Info.SmallestKeptEigenvalue();
         std::cout << '\n';
      }

      MatrixOperator C;
      SL.ConstructMatrices(SL.begin(), Cutoff, *I, C, R);
      R = prod(C, R);

      I=J;
      ++J;
      ++BondIter;
      ++CurrentBond;
   }
   CHECK(BondIter == BondOperators.end());

   *I = R;  // set the final matrix
}

void SweepLeftEvolve(LinearWavefunction& Psi, std::list<SimpleOperator> const& BondOperators,
                     StatesInfo const& SInfo, bool ShowInfo)
{
   LinearWavefunction::iterator J = Psi.end();
   LinearWavefunction::iterator I = J; --I;
   std::list<SimpleOperator>::const_iterator BondIter = BondOperators.end();
   MPStateComponent R = *I;

   int CurrentBond = Psi.size() - 1;
   while (I != Psi.begin())
   {
      CHECK(BondIter != BondOperators.begin());

      J = I;
      --I;
      --BondIter;
      --CurrentBond;

      MPStateComponent A = local_prod(*BondIter, local_tensor_prod(*I, R));
      AMatSVD SL(A, Tensor::ProductBasis<BasisList, BasisList>(I->LocalBasis(), R.LocalBasis()));
      TruncationInfo Info;
      AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(),
                                                                  SInfo, Info);
      if (ShowInfo)
      {
         std::cout << "Bond=(" << (CurrentBond+1) << ',' << (CurrentBond+2) << ")";
         if (ShowStates) std::cout << " states=" << Info.KeptStates();
         if (ShowEntropy) std::cout << " entropy=" << Info.KeptEntropy();
         if (ShowTruncation) std::cout << " trunc=" << Info.TruncationError();
         if (ShowSmallestKept) std::cout << " skeep=" << Info.SmallestKeptEigenvalue();
         std::cout << '\n';
      }

      MatrixOperator C;
      SL.ConstructMatrices(SL.begin(), Cutoff, R, C, *J);
      R = prod(R, C);
   }
   CHECK(BondIter == BondOperators.begin());

   *I = R;  // set the final matrix
}
