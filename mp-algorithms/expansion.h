// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/expansion.h
//
// Copyright (C) 2004-2024 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EXPANSION_H)
#define MPTOOLKIT_MP_ALGORITHMS_EXPANSION_H

#include "mps/state_component.h"
#include "mpo/operator_component.h"
#include "mps/density.h"
#include "mps/truncation.h"
#include "common/namedenum.h"

// a NamedEnumeration for the possible pre-expansion algorithms
struct PreExpansionTraits
{
   enum Enum { SVD, RSVD, RangeFinding, Random, NoExpansion };
   static constexpr std::array<char const*,5> Names = { "fullsvd", "rsvd", "range", "random", "none" };
   static constexpr Enum Default = NoExpansion;
   static constexpr char const* StaticName = "pre-expansion algorithm";
};

using PreExpansionAlgorithm = NamedEnumeration<PreExpansionTraits>;

// a NamedEnumeration for the possible post-expansion algorithms
struct PostExpansionTraits
{
   enum Enum { SVD, RSVD, RangeFinding, Random, Mixing, NoExpansion };
   static constexpr std::array<char const*, 6> Names = { "fullsvd", "rsvd", "range", "random", "mixing", "none"};
   static constexpr Enum Default = RSVD;
   static constexpr char const* StaticName = "post-expansion algorithm";
};

using PostExpansionAlgorithm = NamedEnumeration<PostExpansionTraits>;

// structure to represent the 'oversampling' of the random range finding algorithm for the randomized SVD.
// In each quantum number sector, if the desired number of states is n,
// then over-sample to min(n+Add, n*Scale)
// Recommended value of Add is 10
// Scale should work OK from 1.0 and up.
// ExtraPerSector is an additional 'ExtraStatesPerSector' that applies at the oversampling phase
struct OversamplingInfo
{
   OversamplingInfo() : Add(0), Scale(1.0), ExtraPerSector(0.0) {}
   OversamplingInfo(int Add_, double Scale_) : Add(Add_), Scale(Scale_), ExtraPerSector(0.0) {}
   OversamplingInfo(int Add_, double Scale_, int Extra_) : Add(Add_), Scale(Scale_), ExtraPerSector(Extra_) {}

   // return the actual number of vectors to use, if we want to get k accurate vectors
   int operator()(int k) const { return std::min(k+Add, int(k*Scale+0.5)); }  // round up

   int Add;
   double Scale;
   int ExtraPerSector;
};

std::ostream& operator<<(std::ostream& out, OversamplingInfo const& x);

// Expand the Basis1 of C.
// On exit, Result' * C' = C (up to truncation!), and C is right-orthogonal
MatrixOperator
TruncateExpandBasis1(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector, TruncationInfo& Info, OversamplingInfo Oversampling);


// Apply subspace expansion / truncation on the right (C.Basis2()).
// On exit, C' * Result' = C (up to truncation!), and C' is left-orthogonal
MatrixOperator
TruncateExpandBasis2(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector, TruncationInfo& Info, OversamplingInfo Oversampling);

// Pre-expand Basis 2 of C, by adding vectors to the right-ortho Basis1() of R.
// Returns the expansion vectors in the basis of R.
// Result'.Basis2() == R.Basis2()
// Result'.Basis1() can be anything that doesn't exhaust the rank of R.Basis1() + Result'.Basis1(),
// eg returning random vectors is OK, they do not need to be in the null space as it is assumed that the
// caller will handle that.
StateComponent
PreExpandBasis2(StateComponent const& C, StateComponent const& R, StateComponent const& LeftHam, OperatorComponent const& HLeft, OperatorComponent const& HRight, StateComponent const& RightHam, PreExpansionAlgorithm Algo, int ExtraStates, int ExtraStatesPerSector, OversamplingInfo Oversampling, bool ProjectTwoSiteTangent);

// Pre-expand Basis 1 of C, by adding vectors to the left-ortho Basis2() of L.
// Returns the expansion vectors in the basis of L.
// Result'.Basis1() == L.Basis1()
// Result'.Basis2() can be anything that doesn't exhaust the rank of L.Basis2() + Result'.Basis2(),
// eg returning random vectors is OK, they do not need to be in the null space as it is assumed that the
// caller will handle that.
StateComponent
PreExpandBasis1(StateComponent const& L, StateComponent const& C, StateComponent const& LeftHam, OperatorComponent const& HLeft, OperatorComponent const& HRight, StateComponent const& RightHam, PreExpansionAlgorithm Algo, int ExtraStates, int ExtraStatesPerSector, OversamplingInfo Oversampling, bool ProjectTwoSiteTangent);

// Apply subspace expansion / truncation on the left (C.Basis1()).
// Returns Lambda matrix (diagonal) and a unitary matrix
// Postcondition: U' Lambda' C' = C (up to truncation!)
std::pair<MatrixOperator, RealDiagonalOperator>
SubspaceExpandBasis1(StateComponent& C, OperatorComponent const& H, StateComponent const& RightHam,
                     double MixFactor, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& LeftHam);

// Apply subspace expansion / truncation on the right (C.Basis2()).
// Returns Lambda matrix (diagonal) and a unitary matrix
// Postcondition: C' Lambda' U' = C (up to truncation!)
std::pair<RealDiagonalOperator, MatrixOperator>
SubspaceExpandBasis2(StateComponent& C, OperatorComponent const& H, StateComponent const& LeftHam,
                     double MixFactor, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& RightHam);

#endif
