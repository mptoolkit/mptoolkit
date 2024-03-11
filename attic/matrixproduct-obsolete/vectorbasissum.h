// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/vectorbasissum.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
//
// Alternate version of tensorsum.h
//
// Instead of constructing a SumBasis, we construct the
// projection operators that map the summed basis to
// the constituent parts.  A similar procedure works
// for constructing the truncation operators.
//
// TODO: currently, the projection operators are
// ordinary MatrixOperator's, we could do better than
// that.

#if !defined(MATRIXOPERATORSUM_H_KJDSHY438757438OH)
#define MATRIXOPERATORSUM_H_KJDSHY438757438OH

#include "mpstate.h" // for MatrixOperator definition

typedef MatrixOperator DiagonalProjection;

// Given two VectorBasis, constructs the direct sum
// and returns two projections that map from the sum basis
// to B1 and B2, respectively.
std::pair<DiagonalProjection, DiagonalProjection>
basis_sum(VectorBasis const& B1, VectorBasis const& B2);

// Given a VactorBasis and an iterator pair
// over a container (the kept states)
// that contains Iter->Subspace and
// Iter->Index components, returns projection operators
// onto the kept states and discarded states respectively.
template <typename FwdIter>
std::pair<DiagonalProjection, DiagonalProjection>
basis_truncation(VectorBasis const& B1, FwdIter KeepFirst, FwdIter KeepLast);

// implementation

template <typename FwdIter>
std::pair<DiagonalProjection, DiagonalProjection>
basis_truncation(VectorBasis const& B1, FwdIter KeepFirst, FwdIter KeepLast)
{
   VectorBasis Keep(B1.GetSymmetryList()), Discard(B1.GetSymmetryList());

   // Set of discarded and kept states
   std::vector<std::set<int> > KeepStates(B1.size());
   std::vector<std::set<int> > DiscardStates(B1.size());
   typedef std::set<int>::const_iterator SetIter;
   for (unsigned k = 0; k < B1.size(); ++k)
      for (int l = 0; l < B1.dim(k); ++l)
         DiscardStates[k].insert(l);

   // Make a pass over the kept states to get the set of kept and discarded states
   for (FwdIter I = KeepFirst; I != KeepLast; ++I)
   {
      KeepStates[I->Subspace].insert(I->Index);
      DiscardStates[I->Subspace].erase(I->Index);
   }

   // Now we construct the kept and discarded basis, and the mapping of
   // the subspaces
   std::vector<int> KeepSubspace(B1.size(), -1);
   std::vector<int> DiscardSubspace(B1.size(), -1);
   for (unsigned k = 0; k < B1.size(); ++k)
   {
      if (!KeepStates[k].empty())
      {
         KeepSubspace[k] = Keep.size();
         Keep.push_back(B1[k], KeepStates[k].size());
      }
      if (!DiscardStates[k].empty())
      {
         DiscardSubspace[k] = Discard.size();
         Discard.push_back(B1[k], DiscardStates[k].size());
      }
   }

   // Now we have the kept and discarded basis
   DiagonalProjection KeepP(Keep, B1);
   DiagonalProjection DiscardP(Discard, B1);

   // Finally, one more pass to construct the truncators
   for (unsigned k = 0; k < B1.size(); ++k)
   {
      if (!KeepStates[k].empty())
         KeepP(KeepSubspace[k], k) =
            LinearAlgebra::Matrix<double>(Keep.dim(KeepSubspace[k]), B1.dim(k), 0.0);
      int i = 0;
      for (SetIter I = KeepStates[k].begin(); I != KeepStates[k].end(); ++I,++i)
      {
         KeepP(KeepSubspace[k], k)(i,*I) = 1.0;
      }
      if (!DiscardStates[k].empty())
         DiscardP(DiscardSubspace[k], k) =
            LinearAlgebra::Matrix<double>(Discard.dim(DiscardSubspace[k]), B1.dim(k), 0.0);
      i = 0;
      for (SetIter I = DiscardStates[k].begin(); I != DiscardStates[k].end(); ++I,++i)
      {
         DiscardP(DiscardSubspace[k], k)(i,*I) = 1.0;
      }
   }

   return std::pair<DiagonalProjection, DiagonalProjection>(KeepP, DiscardP);
}

#endif
