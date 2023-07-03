// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/regularize.cpp
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

#include "regularize.h"

namespace Tensor
{

Regularizer::Regularizer(VectorBasis const& b)
{
   if (is_regular_basis(b))
   {
      Trivial = true;
      RegularBasis = b;
      return;
   }

   Trivial = false;
   IrregularBasis = b;
   // Iterate through b and determine the total dimension of
   // each quantum number space, and also map the subspaces of b
   // onto a range of the total space.
   std::map<QuantumNumbers::QuantumNumber, int> SizeMap;
   BasisMappingRange.reserve(b.size());

   for (std::size_t i = 0; i < b.size(); ++i)
   {
      int Sz = SizeMap[b[i]];
      BasisMappingRange.push_back(LinearAlgebra::range(Sz, Sz+b.dim(i)));
      SizeMap[b[i]] += b.dim(i);
   }

   RegularBasis = VectorBasis(SizeMap.begin(), SizeMap.end());

   // <Map the quantum numbers in the basis to the index in the regular basis
   std::map<QuantumNumbers::QuantumNumber, int> IndexOfQ;
   {
      int i = 0;
      for (auto I = SizeMap.begin(); I != SizeMap.end(); ++I)
         IndexOfQ[I->first] = i++;
   }
   BasisMappingIndex.reserve(b.size());
   for (std::size_t i = 0; i < b.size(); ++i)
   {
      BasisMappingIndex.push_back(IndexOfQ[b[i]]);
   }
}

int Regularizer::IndexOf(int i) const
{
   return Trivial ? 1 : BasisMappingIndex[i];
}

LinearAlgebra::Range Regularizer::RangeOf(int i) const
{
   return Trivial ? LinearAlgebra::Range(0, RegularBasis.dim(i)) : BasisMappingRange[i];
}

bool is_regular_basis(VectorBasis const& b)
{
   std::set<QuantumNumber> Used;
   for (std::size_t i = 0; i < b.size(); ++i)
   {
      if (Used.find(b[i]) != Used.end())
         return false;

      Used.insert(b[i]);
   }
   return true;
}

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList>
Regularize(BasisList const& b)
{
   typedef IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList> ResultType;
   // Iterate through b and determine the total dimension of
   // each quantum number space, and also map the subspaces of b
   // onto an index of the total space.
   typedef std::map<QuantumNumbers::QuantumNumber, int> SizeMapType;
   SizeMapType SizeMap;
   std::vector<int> IndexOfSubspace;      // the multiplicity of the quantum number in the basis
   IndexOfSubspace.reserve(b.size());

   for (std::size_t i = 0; i < b.size(); ++i)
   {
      IndexOfSubspace.push_back(SizeMap[b[i]]++);
   }

   // Now make another quantum number map, this time enumerating the quantum numbers in order
   SizeMapType IndexOfQ;
   {
      int i = 0;
      for (SizeMapType::const_iterator I = SizeMap.begin(); I != SizeMap.end(); ++I)
         IndexOfQ[I->first] = i++;
   }

   // Make the final basis
   VectorBasis RegularBasis(SizeMap.begin(), SizeMap.end());

   // and the transform operator
   ResultType Result(RegularBasis, b, QuantumNumbers::QuantumNumber(b.GetSymmetryList()));
   for (std::size_t i = 0; i < b.size(); ++i)
   {
      int Dest = IndexOfQ[b[i]];
      if (size1(Result(Dest, i)) == 0)
         Result(Dest, i) = LinearAlgebra::Matrix<double>(RegularBasis.dim(Dest),1, 0.0);

      Result(Dest,i)(IndexOfSubspace[i], 0) = 1.0;
   }

   return Result;
}

IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList>
SplitBasis(VectorBasis const& b)
{
   BasisList ResultBasis(b.GetSymmetryList());
   for (unsigned i = 0; i < b.size(); ++i)
   {
      for (int j = 0; j < b.dim(i); ++j)
      {
         ResultBasis.push_back(b[i]);
      }
   }

   IrredTensor<LinearAlgebra::Matrix<double>, VectorBasis, BasisList> Result(b, ResultBasis, QuantumNumber(b.GetSymmetryList()));
   int Index = 0;
   for (unsigned i = 0; i < b.size(); ++i)
   {
      for (int j = 0; j < b.dim(i); ++j)
      {
         Result(i, Index) = LinearAlgebra::Matrix<double>(b.dim(i), 1, 0.0);
         Result(i, Index)(j, 0) = 1.0;
         ++Index;
      }
   }
   return Result;
}

} // namespace Tensor
