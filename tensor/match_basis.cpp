// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/match_basis.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "match_basis.h"

namespace Tensor
{

   typedef std::list<std::pair<int, int> > SubIndexPair;
typedef std::map<QuantumNumbers::QuantumNumber, SubIndexPair> SubPairMap;

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
MatchBasis(VectorBasis const& B1, VectorBasis const& B2)
{
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
      Result(B1, B2);

   std::set<QuantumNumbers::QuantumNumber> UsedQN;

   // loop over basis 1.
   SubPairMap Count1;
   for (unsigned i = 0; i < B1.size(); ++i)
   {
      UsedQN.insert(B1[i]);
      for (int j = 0; j < B1.dim(i); ++j)
      {
         Count1[B1[i]].push_back(std::make_pair(i,j));
      }
   }

   // loop over basis 2
   SubPairMap Count2;
   for (unsigned i = 0; i < B2.size(); ++i)
   {
      UsedQN.insert(B2[i]);
      for (int j = 0; j < B2.dim(i); ++j)
      {
         Count2[B2[i]].push_back(std::make_pair(i,j));
      }
   }

   // go through the used quantum numbers and fill the result operator
   for (std::set<QuantumNumbers::QuantumNumber>::const_iterator I = UsedQN.begin();
        I != UsedQN.end(); ++I)
   {
      if (Count1[*I].empty() || Count2[*I].empty()) continue;

      SubIndexPair::const_iterator I1 = Count1[*I].begin();
      SubIndexPair::const_iterator I2 = Count2[*I].begin();

      while (I1 != Count1[*I].end() && I2 != Count2[*I].end())
      {
         if (size1(Result(I1->first, I2->first)) == 0)
            Result(I1->first, I2->first) =
               LinearAlgebra::Matrix<double>(B1.dim(I1->first), B2.dim(I2->first), 0.0);

         Result(I1->first, I2->first)(I1->second, I2->second) = 1.0;

         ++I1;
         ++I2;
      }
   }

   return Result;
}

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
MatchBasisReverse(VectorBasis const& B1, VectorBasis const& B2)
{
   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, VectorBasis>
      Result(B1, B2);

   std::set<QuantumNumbers::QuantumNumber> UsedQN;

   // loop over basis 1.  add pairs to the front, to reverse the sort order
   SubPairMap Count1;
   for (unsigned i = 0; i < B1.size(); ++i)
   {
      UsedQN.insert(B1[i]);
      for (int j = 0; j < B1.dim(i); ++j)
      {
         Count1[B1[i]].push_front(std::make_pair(i,j));
      }
   }

   // loop over basis 2
   SubPairMap Count2;
   for (unsigned i = 0; i < B2.size(); ++i)
   {
      UsedQN.insert(B2[i]);
      for (int j = 0; j < B2.dim(i); ++j)
      {
         Count2[B2[i]].push_front(std::make_pair(i,j));
      }
   }

   // go through the used quantum numbers and fill the result operator
   for (std::set<QuantumNumbers::QuantumNumber>::const_iterator I = UsedQN.begin();
        I != UsedQN.end(); ++I)
   {
      if (Count1[*I].empty() || Count2[*I].empty()) continue;

      SubIndexPair::const_iterator I1 = Count1[*I].begin();
      SubIndexPair::const_iterator I2 = Count2[*I].begin();

      while (I1 != Count1[*I].end() && I2 != Count2[*I].end())
      {
         if (size1(Result(I1->first, I2->first)) == 0)
            Result(I1->first, I2->first) =
               LinearAlgebra::Matrix<double>(B1.dim(I1->first), B2.dim(I2->first), 0.0);

         Result(I1->first, I2->first)(I1->second, I2->second) = 1.0;

         ++I1;
         ++I2;
      }
   }

   return Result;
}

} // namespace Tensor
