// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/dmrg.cpp
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
// $Id$

#include "dmrg.h"
#include "lanczos-ortho.h"
#include "davidson.h"
#include "arnoldi.h"
#include "mps/density.h"
#include "mps/truncation.h"
#include "pstream/optional.h"
#include <boost/optional.hpp>
#include <boost/none.hpp>
#include "linearalgebra/matrix_utility.h"
#include "common/statistics.h"
#include "common/environment.h"
#include "tensor/tensor_eigen.h"
#include <cctype>

double const PrefactorEpsilon = 1e-16;

using MessageLogger::Logger;
using MessageLogger::msg_log;

// These are not needed in C++17 (?), but C++14 requires the initializer
constexpr std::array<char const*,4> PreExpansionTraits::Names;
constexpr std::array<char const*, 5> PostExpansionTraits::Names;

// Helper function to construct a projector that removes the first and last components from a BasisList
SimpleOperator ProjectRemoveFirstLast(BasisList const& B)
{
   SimpleOperator Result(BasisList(B.GetSymmetryList(), B.begin()+1, B.end()-1), B);
   for (int i = 0; i < Result.Basis1().size(); ++i)
   {
      Result(i,i+1) = 1.0;
   }
   return Result;
}

// Helper function to construct a diagonal weight matrix from the frobenius norm of an E/F matrix
SimpleOperator ConstructWeights(StateComponent const& s)
{
   std::vector<double> w;
   w.reserve(s.LocalBasis().size());
   for (auto const& x : s)
      w.push_back(norm_frob(x));

   // If any weights are zero, then manually adjust it to be equal to the smallest non-zero weight.
   // If all of the weights are zero, then set them all to be equal.
   double Min = 0.0;
   for (auto m : w)
   {
      if (m > 0.0 && (m < Min || Min == 0.0))
         Min = m;
   }
   // If all of the weights are identically zero, then set them to be the same
   if (Min == 0.0)
      Min = 1.0;
   for (auto& m : w)
   {
      m = std::max(m, Min);
   }
   SimpleOperator Result(s.LocalBasis());
   for (int i = 0; i < w.size(); ++i)
   {
      Result(i,i) = w[i];
   }
   return Result;
}

// Helper function to distribute ExtraStates among the available states, with relative weights given by the Weights array.
// To do this, we start from the total weight of the available sectors, and stretch that interval w[j] to be length ExtraStates.
// (we can cap ExtraStates at the total number of available states in sectors with non-zero weight, so ExtraStates is actually
// achievable. We can also remove sectors that have zero weight from the Avail array.)
// Then the number of states we keep of the j'th sector, k[j] is k[j] = round(sum(k<=j) w[k] - k[j-1]) (with k[-1] = 0).
// If any of the k[j] > Avail[j], then pin k[j] = Avail[j], remove sector j from the Avail array, and repeat the distribution.
// We avoid creating subspaces that have zero dimension. In principle this is harmless, although inefficient, but
// MKL doesn't behave properly when multiplying 0-size matrices.
// CapExtraStatesPerSector is a cutoff in case there is a large number of different quantum number sectors available.
// Set it to zero for no limit, but if it is a positive number then it limits the total number of states added by
// ExtraStatesPerSector (in which case they are chosen at random).
std::map<QuantumNumbers::QuantumNumber, int>
DistributeStates(std::map<QuantumNumbers::QuantumNumber, int> Weights, std::map<QuantumNumbers::QuantumNumber, int> const& Avail, int ExtraStates, int ExtraStatesPerSector, int CapExtraStatesPerSector = 0)
{
   std::map<QuantumNumbers::QuantumNumber, int> Result;

   // Determine the total count of available states in sectors that have non-zero weight, and their total weight
   int TotalWeight = 0;
   int TotalAvail = 0;
   for (auto const& a : Avail)
   {
      if (Weights[a.first] > 0)
      {
         TotalWeight += Weights[a.first];
         TotalAvail += a.second;
      }
   }
   ExtraStates = std::min(ExtraStates, TotalAvail); // cap ExtraStates, if it is more than the states available

   // Fill the Result array with the interval length of the Weights, stretched to length ExtraStates
   bool Valid = false;
   while (!Valid)
   {
      Valid = true;
      int StatesKeptSoFar = 0;
      int WeightSum = 0;
      int TotalWeightThisRound = TotalWeight;
      int StatesWantedThisRound = ExtraStates;
      for (auto const& a : Avail)
      {
         if (Weights[a.first] > 0)
         {
            WeightSum += Weights[a.first];
            int NumToAdd = int(std::round((double(WeightSum) / TotalWeightThisRound) * StatesWantedThisRound)) - StatesKeptSoFar;
            if (NumToAdd > 0)
            {
               Result[a.first] = NumToAdd;
               StatesKeptSoFar += Result[a.first];

               // if we've exceed the possible allocation, then peg the number of states at the maximum and
               // remove this sector from consideration in the next pass
               if (NumToAdd > a.second)
               {
                  Result[a.first] = a.second;
                  ExtraStates -= a.second;
                  TotalWeight -= Weights[a.first];
                  Weights[a.first] = 0;
                  Valid = false;
               }
            }
         }
      }
      DEBUG_CHECK_EQUAL(StatesKeptSoFar, StatesWantedThisRound);
   }

   if (ExtraStatesPerSector > 0)
   {
      // If we don't have a cap on the number of states per sector to add, then just go through the list and add them
      if (CapExtraStatesPerSector == 0)
      {
         // Final pass: add ExtraStatesPerSector to states that don't already have them
         for (auto a : Avail)
         {
            if (Result[a.first] < ExtraStatesPerSector && Result[a.first] < a.second)
               Result[a.first] = std::min(a.second, ExtraStatesPerSector);
         }
      }
      else
      {
         // If we have a cap on the number, we need to select states to keep at random.
         // Firstly assemble a list of all available sectors, including multiple times if ExtraStatesPerSector > 1
         std::vector<QuantumNumbers::QuantumNumber> ToAdd;
         for (auto a : Avail)
         {
            int AddedThisSector = Result.count(a.first) == 0 ? 0 : Result[a.first];
            if (AddedThisSector < ExtraStatesPerSector && AddedThisSector < a.second)
            {
               for (int i = 0; i < std::min(a.second, ExtraStatesPerSector) - AddedThisSector; ++i)
               {
                  ToAdd.push_back(a.first);
               }
            }
         }
         // randomize the order of the sectors, but don't bother if we're going to add all of the states anyway
         if (ToAdd.size() > CapExtraStatesPerSector)
         {
            randutil::random_stream stream;
            stream.seed();
            std::shuffle(ToAdd.begin(), ToAdd.end(), stream.u_rand);
         }
         // Now add states up to the cap
         for (auto a = ToAdd.begin(); a != ToAdd.end() && CapExtraStatesPerSector > 0; ++a)
         {
            ++Result[*a];
            --CapExtraStatesPerSector;
         }
      }
   }
   return Result;
}

// Get a basis for use in the randomized SVD, selecting ExtraStates*RangeFindingOverhead states out of the available states.
// The quantum number sectors for the singular vectors are chosen according to the relative number of states in
// KeptStateDimension, with at least ExtraStatesPerSector states in each available quantum number sector.
// RangeFindingOverhead is a number >= 1 which indicates how large the embedding space should be for the randomized SVD.
// Typically RangeFindingOverhead = 2.
// Typically KeptStateDimension will be initialized as DimensionPerSector() of some basis of kept states.
// On entry: X is a p x q matrix (typically = p m*d, q = w*m)
// On exit: Result' is a VectorBasis of size ExtraStates + additional states per sector
VectorBasis
MakeExpansionBasis(QuantumNumbers::SymmetryList const& SL, std::map<QuantumNumbers::QuantumNumber, int> const& NumAvailablePerSector, std::map<QuantumNumber, int> const& KeptStateDimension, int ExtraStates, int ExtraStatesPerSector, double RangeFindingOverhead)
{
   std::map<QuantumNumbers::QuantumNumber, int> Weights;
   for (auto n : NumAvailablePerSector)
   {
      // Set the relative weight to the number of kept states in this sector.
      // Clamp this at a minimum of 1 state, since a state only appears in X if it has non-zero contribution
      // with respect to the density matrix mixing.
      if (n.second > 0)
      {
         if (KeptStateDimension.count(n.first))
            Weights[n.first] = std::max(KeptStateDimension.at(n.first), 1);
         else
            Weights[n.first] = 1;
      }
   }
   auto NumExtraPerSector = DistributeStates(Weights, NumAvailablePerSector, int(ExtraStates*RangeFindingOverhead+0.5), int(ExtraStatesPerSector*RangeFindingOverhead+0.5));

   return VectorBasis(SL, NumExtraPerSector.begin(), NumExtraPerSector.end());
}

// Helper function to use the random range-finding algorithm to get the important left-singular vectors of X.
// The candidate basis is obtained using the relative sizes contained in KeptStateDimension.
MatrixOperator RangeFindingBasis1(MatrixOperator X, std::map<QuantumNumber, int> const& KeptStateDimension, int ExtraStates, int ExtraStatesPerSector, double RangeFindingOverhead)
{
   VectorBasis RangeBasis = MakeExpansionBasis(X.GetSymmetryList(), RankPerSector(X), KeptStateDimension, ExtraStates, ExtraStatesPerSector, RangeFindingOverhead);
   MatrixOperator A = X * MakeRandomMatrixOperator(X.Basis2(), RangeBasis);
   auto QR = QR_Factorize(A);
   X = herm(QR.first) * X;

   CMatSVD ExpandDM(X, CMatSVD::Left);

   auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);

   return QR.first * ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
}

// Expand the Basis1 of C.
// On exit, Result' * C' = C (up to truncation!), and C is right-orthogonal
MatrixOperator
TruncateExpandBasis1(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector, TruncationInfo& Info, double RangeFindingOverhead)
{
   // Reshape C to m x dm
   MatrixOperator CMat = ReshapeBasis2(C);

   // R is the right-orthogonal A-matrix of states to keep.  The algorithms below obtain R in different ways.
   StateComponent R;

   // If we don't need to increase the basis size then we can do a 'fast path' that avoids constructing the null space.
   // We can do that if there is no basis expansion and (either the number of states to keep is not larger than the
   // current number of states, OR the StatesInfo specifies that we don't keep zero eigenvalues).
   // Question of what to do if ExtraStates == 0, but ExtraStatesPerSector > 0.  We interpret this to mean that
   // we only want to keep any extra states all if ExtraStates > 0, so we do nothing in this case.
   if (ExtraStates <= 0 && !(States.KeepZeroEigenvalues() && States.MaxStates > C.Basis2().total_dimension()))
   {
      Algo = PostExpansionAlgorithm::NoExpansion;
   }

   if (Algo == PostExpansionAlgorithm::NoExpansion || Algo == PostExpansionAlgorithm::FastRangeFinding)
   {
      // We can do a thin SVD here.  In the no expansion case we just do the basic truncation and we are done.
      // In the fastrangefinding algorithm we implicitly use the full space of m*d states, but project out
      // the kept states.  The arrangement of contractions means that this is faster than constructing the
      // null space explicitly.
      CMatSVD DM(CMat);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
      R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());

      if (Algo == PostExpansionAlgorithm::NoExpansion)
      {
         Info.ExtraStates_ = 0;
      }
      else if (Algo == PostExpansionAlgorithm::FastRangeFinding)
      {
         SimpleOperator P = ProjectRemoveFirstLast(LeftHam.LocalBasis());
         SimpleOperator W = ConstructWeights(local_prod(P, LeftHam));

         // Construct the embedding matrix from the w*m dimensional basis to k dimensions. Firstly get the w*m basis
         VectorBasis RBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(W.Basis1(), C.Basis1()).Basis()).Basis();

         // Get the basis of candidate states based on the ranks of the kept states
         VectorBasis EBasis = MakeExpansionBasis(VKeep.GetSymmetryList(), RankPerSector(VKeep.Basis2(), RBasis), DimensionPerSector(VKeep.Basis1()), ExtraStates, ExtraStatesPerSector, RangeFindingOverhead);

         // Construct the gaussian random embedding matrix
         MatrixOperator M = MakeRandomMatrixOperator(EBasis, RBasis);
         StateComponent Ms = ReshapeFromBasis2(M, W.Basis1(), C.Basis1());

         // The embedding matrix replaces the E matrix elements; the contraction is otherwise the same as a Hamiltonian-vector
         // multiply of the wavefunction.
         StateComponent RNull = operator_prod_inner(W*P*H, Ms, C, herm(RightHam));

         // project out the kept states
         RNull = RNull - scalar_prod(RNull, herm(R)) * R;

         // Range-finding algorithm via LQ decomposition.  We can throw away 'L' here.  It is important that we
         // do the QR *after* projecting out R.  Otherwise we would be doing rangefinding onto the leading
         // singular values including the kept states.
         OrthogonalizeBasis1_LQ(RNull);

         // From here it is the same as the SVD algorithm: construct the mixing term in the projected RNull space
         // and find the expansion vectors from the SVD. The only difference is that RNull will contain some numerical
         // noise from projecting out the kept states, so we need to do an additional explicit orthogonalization.
         // I think we could do that orthogonalization at any time, but slightly faster to do it after projecting
         // down to the final expansion size.
         StateComponent D = contract_from_right(herm(W*P*H), RNull, RightHam, herm(C));
         MatrixOperator X = ReshapeBasis2(D);

         CMatSVD ExpandDM(X, CMatSVD::Left);
         auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
         MatrixOperator QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());

         RNull = herm(QExpand) * RNull;
         Info.ExtraStates_ = RNull.Basis1().total_dimension();

         R = RegularizeBasis1(tensor_col_sum(R, RNull));
         OrthogonalizeBasis1_LQ(R);  // additional orthogonalization step, ensure vectors are all orthogonal
      }
   }
   else if (Algo == PostExpansionAlgorithm::SVD || Algo == PostExpansionAlgorithm::RangeFinding || Algo == PostExpansionAlgorithm::Random)
   {
      // The SVD, rangefinding, and randokm algorithms differ only by rangefinding using a randomized SVD.
      // We need to calculate the null space, so we want all of the right singular vectors
      CMatSVD DM(CMat, CMatSVD::Right);

      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
      MatrixOperator VDiscard = DM.ConstructRightVectors(DMPivot, DM.end());
      R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());         // kept states in right ortho form
      StateComponent RNull = ReshapeFromBasis2(VDiscard, C.LocalBasis(), C.Basis2());  // discarded states in right ortho form

      // Get QExpand, which is the projection from (discarded_states, states_to_expand)
      MatrixOperator QExpand;
      if (Algo == PostExpansionAlgorithm::Random)
      {
         // construct QExpand as Haar random vectors.  Construct the quantum number sectors based on the
         // distribution of kept states. One complication here is how to handle ExtraStatesPerSector. Since we have
         // no idea whether quantum number sectors will be important or not, we limit the number of possible additional
         // sectors to be no more than x times the size of the kept sectors.
         auto NumAvailablePerSector = DimensionPerSector(VDiscard.Basis1());
         auto Weights = DimensionPerSector(VKeep.Basis1());
         auto NumExtraPerSector = DistributeStates(Weights, NumAvailablePerSector, ExtraStates, ExtraStatesPerSector, int(VKeep.Basis1().size() * 1));

         // Make a basis from the distribution of states and get a random matrix between that basis and the discarded states
         VectorBasis ExpansionBasis(C.GetSymmetryList(), NumExtraPerSector.begin(), NumExtraPerSector.end());
         MatrixOperator M = MakeRandomMatrixOperator(VDiscard.Basis1(), ExpansionBasis);

         // Project the discarded basis onto the random vectors via the QR decomposition
         auto QR = QR_Factorize(M);
         QExpand = std::move(QR.first);
      }
      else
      {
         // For SVD and rangefinding, we construct the mixing term from the 3S algorithm
         SimpleOperator P = ProjectRemoveFirstLast(LeftHam.LocalBasis());
         SimpleOperator W = ConstructWeights(local_prod(P, LeftHam));
         MatrixOperator X = ReshapeBasis2(contract_from_right(herm(W*P*H), RNull, RightHam, herm(C)));
         // X is now (discarded basis) x (wm) matrix

         if (Algo == PostExpansionAlgorithm::SVD)
         {
            CMatSVD ExpandDM(X, CMatSVD::Left);
            auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
            QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
         }
         else if (Algo == PostExpansionAlgorithm::RangeFinding)
         {
            QExpand = RangeFindingBasis1(X, DimensionPerSector(VKeep.Basis1()), ExtraStates, ExtraStatesPerSector, RangeFindingOverhead);
         }
      }

      RNull = herm(QExpand) * RNull;
      Info.ExtraStates_ = RNull.Basis1().total_dimension();

      R = RegularizeBasis1(tensor_col_sum(R, RNull));
   }
   else
   {
      PANIC("Unsupported expansion algorithm");
   }

   // Construct the new Lambda matrix
   MatrixOperator Lambda = scalar_prod(C, herm(R));
   C = std::move(R);
   CHECK_EQUAL(LeftHam.Basis1(), Lambda.Basis1());
   CHECK_EQUAL(Lambda.Basis2(), C.Basis1());
   CHECK_EQUAL(C.Basis2(), RightHam.Basis1());
   return Lambda;
}

// Apply subspace expansion / truncation on the right (C.Basis2()).
// On exit, C' * Result' = C (up to truncation!), and C' is left-orthogonal
MatrixOperator
TruncateExpandBasis2(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector, TruncationInfo& Info, double RangeFindingOverhead)
{
   // Reshape C to dm x m
   MatrixOperator CMat = ReshapeBasis1(C);

   // L is the left-orthogonal A-matrix of states to keep.  The algorithms below obtain L in different ways.
   StateComponent L;

   // If we don't need to increase the basis size then we can do a 'fast path' that avoids constructing the null space.
   // We can do that if there is no basis expansion and (either the number of states to keep is not larger than the
   // current number of states, OR the StatesInfo specifies that we don't keep zero eigenvalues).
   if (ExtraStates <= 0 && !(States.KeepZeroEigenvalues() && States.MaxStates > C.Basis2().total_dimension()))
   {
      Algo = PostExpansionAlgorithm::NoExpansion;
   }

   if (Algo == PostExpansionAlgorithm::NoExpansion || Algo == PostExpansionAlgorithm::FastRangeFinding)
   {
      // We can do a thin SVD here.  In the no expansion case we just do the basic truncation and we are done.
      // In the fastrangefinding algorithm we implicitly use the full space of m*d states, but project out
      // the kept states.  The arrangement of contractions means that this is faster than constructing the
      // null space explicitly.
      CMatSVD DM(CMat);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator UKeep = DM.ConstructLeftVectors(DM.begin(), DMPivot);
      L = ReshapeFromBasis1(UKeep, C.LocalBasis(), C.Basis1());

      if (Algo == PostExpansionAlgorithm::FastRangeFinding)
      {
         SimpleOperator P = ProjectRemoveFirstLast(RightHam.LocalBasis());
         SimpleOperator W = ConstructWeights(local_prod(P, RightHam));

         // Construct the embedding matrix from the w*m dimensional basis to k dimensions. Firstly get the w*m basis
         VectorBasis RBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(W.Basis2(), C.Basis2()).Basis()).Basis();

         // Get the basis of candidate states based on the ranks of the kept states
         VectorBasis EBasis = MakeExpansionBasis(UKeep.GetSymmetryList(), RankPerSector(UKeep.Basis1(), RBasis), DimensionPerSector(UKeep.Basis2()), ExtraStates, ExtraStatesPerSector, RangeFindingOverhead);

         // Construct the gaussian random embedding matrix
         MatrixOperator M = MakeRandomMatrixOperator(EBasis, RBasis);
         StateComponent Ms = ReshapeFromBasis2(M, W.Basis2(), C.Basis2());

         // The embedding matrix replaces the F matrix elements; the contraction is otherwise the same as a Hamiltonian-vector
         // multiply of the wavefunction.
         StateComponent LNull = operator_prod_inner(H*herm(P)*W, LeftHam, C, herm(Ms));

         // Project out the kept states
         LNull = LNull - L * scalar_prod(herm(L), LNull);

         // Range-finding algorithm via QR decomposition.  We can throw away 'R' here.  It is important that we
         // do the QR *after* projecting out L.  Otherwise we would be doing rangefinding onto the leading
         // singular values including the kept states.
         // LNull is now (dm) x 2k orthogonal, and is the sample space for the expansion vectors
         OrthogonalizeBasis2_QR(LNull); // throw away the 'R', we don't need it

         // From here it is the same as the SVD algorithm: construct the mixing term in the projected LNull space
         // and find the expansion vectors from the SVD. The only difference is that LNull will contain some numerical
         // noise from projecting out the kept states, so we need to do an additional explicit orthogonalization.
         // I think we could do that orthogonalization at any time, but slightly faster to do it after projecting
         // down to the final expansion size.
         StateComponent D = contract_from_left(H*herm(P)*W, herm(LNull), LeftHam, C);
         MatrixOperator X = ReshapeBasis2(D);

         CMatSVD ExpandDM(X, CMatSVD::Left);
         auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
         MatrixOperator QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
         // UExpand is (discarded_states, states_to_expand)

         LNull = LNull * QExpand;
         Info.ExtraStates_ = LNull.Basis2().total_dimension();

         // LNull might not be exactly orthogonal
         L = RegularizeBasis2(tensor_row_sum(L, LNull));
         OrthogonalizeBasis2_QR(L);
      }
   }
   else if (Algo == PostExpansionAlgorithm::SVD || Algo == PostExpansionAlgorithm::RangeFinding || Algo == PostExpansionAlgorithm::Random)
   {
      // The SVD, rangefinding, and randokm algorithms differ only by rangefinding using a randomized SVD.
      // We need to calculate the null space, so we want all of the left singular vectors
      CMatSVD DM(CMat, CMatSVD::Left);

      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator UKeep = DM.ConstructLeftVectors(DM.begin(), DMPivot);
      MatrixOperator UDiscard = DM.ConstructLeftVectors(DMPivot, DM.end());
      L = ReshapeFromBasis1(UKeep, C.LocalBasis(), C.Basis1());
      StateComponent LNull = ReshapeFromBasis1(UDiscard, C.LocalBasis(), C.Basis1());

      // Get QExpand, which is the projection from (discarded_states, states_to_expand)
      MatrixOperator QExpand;
      if (Algo == PostExpansionAlgorithm::Random)
      {
         // construct QExpand as Haar random vectors.  Construct the quantum number sectors based on the
         // distribution of kept states. One complication here is how to handle ExtraStatesPerSector. Since we have
         // no idea whether quantum number sectors will be important or not, we limit the number of possible additional
         // sectors to be no more than x times the size of the kept sectors.
         auto NumAvailablePerSector = DimensionPerSector(UDiscard.Basis2());
         auto Weights = DimensionPerSector(UKeep.Basis2());
         auto NumExtraPerSector = DistributeStates(Weights, NumAvailablePerSector, ExtraStates, ExtraStatesPerSector, int(UKeep.Basis2().size() * 1));

         // Make a basis from the distribution of states and get a random matrix between that basis and the discarded states
         VectorBasis ExpansionBasis(L.GetSymmetryList(), NumExtraPerSector.begin(), NumExtraPerSector.end());
         MatrixOperator M = MakeRandomMatrixOperator(UDiscard.Basis2(), ExpansionBasis);

         // Project the discarded basis onto the random vectors via the QR decomposition
         auto QR = QR_Factorize(M);
         QExpand = std::move(QR.first);
      }
      else
      {
         // For SVD and rangefinding, we construct the mixing term from the 3S algorithm
         SimpleOperator P = ProjectRemoveFirstLast(RightHam.LocalBasis());
         SimpleOperator W = ConstructWeights(local_prod(P, RightHam));
         MatrixOperator X = ReshapeBasis2(contract_from_left(H*herm(P)*W, herm(LNull), LeftHam, C));
         // X is now (discarded basis) x (wm) matrix

         if (Algo == PostExpansionAlgorithm::SVD)
         {
            // Get the kept and discarded states
            CMatSVD ExpandDM(X, CMatSVD::Left);
            auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
            QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
         }
         else if (Algo == PostExpansionAlgorithm::RangeFinding)
         {
            QExpand = RangeFindingBasis1(X, DimensionPerSector(UKeep.Basis2()), ExtraStates, ExtraStatesPerSector, RangeFindingOverhead);
         }
      }

      LNull = LNull * QExpand;
      Info.ExtraStates_ = LNull.Basis2().total_dimension();

      L = RegularizeBasis2(tensor_row_sum(L, LNull));
   }
   else
   {
      PANIC("Unsupported expansion algorithm");
   }

   MatrixOperator Lambda = scalar_prod(herm(L), C);
   C = std::move(L);
   CHECK_EQUAL(LeftHam.Basis1(), C.Basis1());
   CHECK_EQUAL(C.Basis2(), Lambda.Basis1());
   CHECK_EQUAL(Lambda.Basis2(), RightHam.Basis1());
   return Lambda;
}

PStream::opstream& operator<<(PStream::opstream& out, DMRG const& d)
{
   return out << d.HamMatrices
              << d.Psi
	      << d.Site
              << d.Hamiltonian
	      << d.LeftStop
	      << d.RightStop

              << d.LastOverlap
              << d.IsPsiConverged
              << d.IsConvergedValid
              << d.KeepList

              << d.TotalSweepNumber
              << d.TotalSweepRecNumber
              << d.TotalNumIterations
              << d.TotalNumMultiplies

              << d.SweepNumIterations
              << d.SweepSumStates
              << d.SweepMaxStates
              << d.SweepNumMultiplies
              << d.SweepEnergy
              << d.SweepTruncation
              << d.SweepEntropy
              << d.SweepStartTime
              << d.SweepTruncatedEnergy
              << d.SweepEnergyError
              << d.SweepLastMixFactor

              << d.IterationNumMultiplies
              << d.IterationNumStates
              << d.IterationEnergy
              << d.IterationTruncation
              << d.IterationEntropy
              << d.IterationEnergyBeforeTrunc
              << d.IterationEnergyVec;
}

PStream::ipstream& operator>>(PStream::ipstream& in, DMRG& d)
{
   return in >> d.HamMatrices
             >> d.Psi
	     >> d.Site
             >> d.Hamiltonian
	     >> d.LeftStop
	     >> d.RightStop

             >> d.LastOverlap
             >> d.IsPsiConverged
             >> d.IsConvergedValid
             >> d.KeepList

             >> d.TotalSweepNumber
             >> d.TotalSweepRecNumber
             >> d.TotalNumIterations
             >> d.TotalNumMultiplies

             >> d.SweepNumIterations
             >> d.SweepSumStates
             >> d.SweepMaxStates
             >> d.SweepNumMultiplies
             >> d.SweepEnergy
             >> d.SweepTruncation
             >> d.SweepEntropy
             >> d.SweepStartTime
             >> d.SweepTruncatedEnergy
             >> d.SweepEnergyError
             >> d.SweepLastMixFactor

             >> d.IterationNumMultiplies
             >> d.IterationNumStates
             >> d.IterationEnergy
             >> d.IterationTruncation
             >> d.IterationEntropy
             >> d.IterationEnergyBeforeTrunc
             >> d.IterationEnergyVec;

   // advance the C iterator
   int n = d.Site;
   d.C = d.Psi.begin();
   d.H = d.Hamiltonian.begin();
   while (n > 0)
   {
      ++d.C;
      ++d.H;
      --n;
   }
}

DMRG::DMRG(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_)
   : Hamiltonian(Ham_),
     NormalizeWavefunction(false),
     IsPsiConverged(false), IsConvergedValid(false),
     MixUseEnvironment(false), UseDGKS(false), Solver_(LocalEigensolver::Solver::Lanczos),
     Verbose(Verbose_),
     RangeFindingOverhead(2.0)
{
   // construct the HamMatrix elements
   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;
   H = Hamiltonian.begin();
   HamMatrices.push_left(Initial_E(Ham_, Psi_.Basis1()));
   for (FiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "site " << (HamMatrices.size_left()) << std::endl;
      HamMatrices.push_left(contract_from_left(*H, herm(*I), HamMatrices.left(), *I));
      Psi.push_back(*I);
      ++H;
   }
   HamMatrices.push_right(Initial_F(Ham_, Psi_.Basis2()));
   // Probably no need to incorporate lambda_r(), this is just a normalization
   Psi.set_back(prod(Psi.get_back(), Psi_.lambda_r()));

   // initialize to the right-most site
   HamMatrices.pop_left();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size()-1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;

   // clear the global statistics
   TotalSweepNumber = 0;
   TotalSweepRecNumber = 0;
   TotalNumIterations = 0;
   TotalNumMultiplies = 0;

  this->debug_check_structure();
}

FiniteWavefunctionLeft
DMRG::Wavefunction() const
{
   return FiniteWavefunctionLeft::Construct(Psi);
}

void DMRG::StartIteration()
{
   IterationNumMultiplies = 0;
   IterationNumStates = 0;
   IterationEnergy = 0;
   IterationTruncation = 0;
   IterationEntropy = 0;
}

void DMRG::EndIteration()
{
   // recalculate the energy, as it will have changed after truncation
   if (IterationTruncation != 0)
      IterationEnergy = this->Energy();
   msg_log(1, "EnergyLog") << TotalSweepNumber << ' '
                           << TotalSweepRecNumber << ' '
                           << Site << ' '
                           << IterationNumStates << ' '
                           << IterationNumMultiplies << ' '
                           << IterationTruncation << ' '
                           << IterationEntropy << ' '
                           << IterationEnergy << '\n';

   ++SweepNumIterations;
   SweepNumMultiplies += IterationNumMultiplies;
   SweepSumStates += IterationNumStates;
   SweepMaxStates = std::max(SweepMaxStates, IterationNumStates);
   SweepEnergy = std::min(SweepEnergy, IterationEnergy.real());
   SweepTruncatedEnergy += IterationEnergyBeforeTrunc.real() - IterationEnergy.real();
   SweepTruncation += IterationTruncation;
   SweepEntropy = std::max(SweepEntropy, IterationEntropy);
   IterationEnergyVec.push_back(IterationEnergy);
}

void DMRG::StartSweep(bool IncrementSweepNumber, double /* Broad_ */)
{
   ++TotalSweepNumber;
   if (IncrementSweepNumber)
      ++TotalSweepRecNumber;

   SweepNumIterations = 0;
   SweepSumStates = 0;
   SweepMaxStates = 0;
   SweepNumMultiplies = 0;
   SweepEnergy = 1E100;
   SweepTruncation = 0;
   SweepEntropy = 0;
   SweepStartTime = ProcControl::GetCumulativeElapsedTime();
   SweepTruncatedEnergy = 0;
   IterationEnergyVec.clear();

   PsiPrevC = *C;

   // Initialize the keep list
   KeepList.clear();

   IsConvergedValid = false;
}

double DMRG::FidelityLoss() const
{
#if 0
   return 2.0 * (1.0 - norm_frob(inner_prod(Psi.Center(), PsiPrevC))
                 / norm_frob(Psi.Center()));
#endif
}

void DMRG::EndSweep()
{
   TotalNumIterations += SweepNumIterations;
   TotalNumMultiplies += SweepNumMultiplies;

   SweepEnergyError = std::sqrt(statistics::variance(IterationEnergyVec.begin(),
                                                     IterationEnergyVec.end())
                                / SweepNumIterations);

   double Overlap = 2.0 * (1.0 - norm_frob(inner_prod(*C, PsiPrevC))
                           / norm_frob(*C));
   double OverlapDifference = 2;
   if (LastOverlap)
   {
      OverlapDifference = std::abs(*LastOverlap - Overlap);
      IsPsiConverged = ((Overlap < ConvergenceOverlapTruncationScale * SweepTruncation)
                        || SweepTruncation < ConvergenceSweepTruncMin)
         && ((OverlapDifference < ConvergenceOverlapDifferenceOverlapScale * Overlap)
             || Overlap < ConvergenceOverlapMin);
      IsConvergedValid = true;
   }

   LastOverlap = Overlap;


   int Converged = this->IsConverged();
   double SweepTime = ProcControl::GetCumulativeElapsedTime() - SweepStartTime;

   msg_log(1, "SweepLog") << TotalSweepNumber << ' '
                          << TotalSweepRecNumber << ' '
                          << SweepNumIterations << ' '
                          << SweepMaxStates << ' '
                          << SweepNumMultiplies << ' '
                          << SweepTime << ' '
                          << std::max(SweepTruncation, 0.0) << ' '
                          << SweepTruncatedEnergy << ' '
                          << SweepEnergyError << ' '
                          << SweepEntropy << ' '
                          << Overlap << ' '
                          << OverlapDifference << ' '
                          << SweepEnergy << ' '
                          << SweepLastMixFactor << ' '
                          << Converged
                          << std::endl;
   msg_log(1, "EnergyLog") << std::flush;
}

std::complex<double>
DMRG::Energy() const
{
   return inner_prod(HamMatrices.right(),
		     contract_from_left(*H, herm(HamMatrices.left()), *C, HamMatrices.right()));
}

std::complex<double>
DMRG::Solve()
{
   Solver_.Solve(*C, HamMatrices.left(), *H, HamMatrices.right());

   // NOTE: *C is typically nearly, but not exactly normalized at this point.
   // There is no point normalizing here, since we need to do it anyway after the truncation.

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   IterationEnergyBeforeTrunc = IterationEnergy;

   return IterationEnergy;
}

bool DMRG::IsConverged() const
{
   return IsConvergedValid && IsPsiConverged;
}

SimpleOperator ProjectMiddle(BasisList const& B)
{
   SimpleOperator Result(B, BasisList(B.GetSymmetryList(), B.begin()+1, B.end()-1));
   int Size = Result.Basis2().size();
   for (int i = 0; i < Size; ++i)
   {
      Result(i+1,i) = 1.0;
   }
      return Result;
}

// Given a two-site section of MPS, given by (L, C) in left and center ortho respectively, expand the basis between them.
// E and F are the Hamiltonian matrices for L.Basis1() and C.Basis2() respectively.
// Returns the number of states that were added to the environment basis.
int
ExpandLeftEnvironment(StateComponent& L, StateComponent& C,
                      StateComponent const& E, StateComponent const& F,
                      OperatorComponent const& HL, OperatorComponent const& HR,
                      int StatesWanted, int ExtraStatesPerSector)
{
   DEBUG_CHECK_EQUAL(L.Basis2(), C.Basis1());
   DEBUG_CHECK_EQUAL(E.Basis1(), L.Basis1());
   DEBUG_CHECK_EQUAL(F.Basis1(), C.Basis2());
   DEBUG_CHECK_EQUAL(L.LocalBasis(), HL.LocalBasis2());
   DEBUG_CHECK_EQUAL(C.LocalBasis(), HR.LocalBasis2());

#if 0
   int ExtraStates = StatesWanted - C.Basis1().total_dimension();
   // Calculate left null space of left site.
   StateComponent LNull = NullSpace2(L);

   // Perform SVD to right-orthogonalize the right site and extract singular value matrix.
   StateComponent R = C;
   MatrixOperator U;
   RealDiagonalOperator LambdaD;
   std::tie(U, LambdaD) = OrthogonalizeBasis1(R); // MPS is L,U,Lambda,R
   MatrixOperator Lambda = U*LambdaD;

   auto CLeft = L*U*LambdaD;

   StateComponent RightHam = contract_from_right(herm(HR), R, F, herm(R));

   StateComponent = DensityMatrixMixingBasis2(LNull, CLeft, LambdaD, E, HL, RightHam);

   // Now the 3S-inspired step: calculate X = new F matrix projected onto the null space.
   // Firstly project out the first and last columns of HLeft.
   // NOTE: if the Hamiltonian is 2-dimensional MPO (i.e. there are no interactions here), then
   // the projector maps onto an empty set, so we need the 3-parameter constructor of the BasisList
   SimpleOperator Projector(HL.Basis2(), BasisList(HLeft.GetSymmetryList(), HL.Basis2().begin()+1, HLeft.Basis2().end()-1));
   for (int i = 0; i < Projector.Basis2().size(); ++i)
      Projector(i+1, i) = 1.0;



   StateComponent X = contract_from_left(HL*Projector, herm(LNull), E, CL*URight*DRight);
   StateComponent FRight = contract_from_right(herm(herm(Projector)*HRight), COrtho, F, herm(C));

   // Multiply each element of X by a prefactor depending on the corresponding element of F.
   for (int i = 0; i < X.size(); ++i)
   {
      // We add an epsilon term to the prefactor to avoid the corner case where
      // the prefactor is zero for all elements, and any possible new states
      // are ignored.
      double Prefactor = norm_frob(FRight[i]) + PrefactorEpsilon;
      X[i] *= Prefactor;
   }

   // Take the truncated SVD of X.
   MatrixOperator XExpand = ReshapeBasis2(X);
   XExpand = MatrixOperator(XExpand.Basis1(), VectorBasis(XExpand.GetSymmetryList()));

   CMatSVD SVD(XExpand, CMatSVD::Left);

   auto StatesToKeep = TruncateExtraStates(SVD.begin(), SVD.end(), ExtraStates, ExtraStatesPerSector, true);

   MatrixOperator U = SVD.ConstructLeftVectors(StatesToKeep.begin(), StatesToKeep.end());

   // Construct new basis.
   SumBasis<VectorBasis> NewBasis(CLeft.Basis2(), U.Basis2());
   // Regularize the new basis.
   Regularizer R(NewBasis);

   // Add the new states to CLeft, and add zeros to C.
   CLeft = RegularizeBasis2(tensor_row_sum(CLeft, prod(NLeft, U), NewBasis), R);

   StateComponent Z = StateComponent(C.LocalBasis(), U.Basis2(), C.Basis2());
   C = RegularizeBasis1(R, tensor_col_sum(C, Z, NewBasis));

   return StatesToKeep.size();
   #endif
}

// helper to subtract one std::Map<T,int> from another
template <typename T>
std::map<T, int>
subtract(std::map<T, int> x, std::map<T, int> const& y)
{
   for (auto const& yi : y)
   {
      if ((x[yi.first] -= yi.second) == 0)
         x.erase(yi.first);
   }
   return x;
}

// Pre-expand Basis 2 of C, by adding vectors to the right-ortho Basis1() of R.
void
PreExpandBasis2(StateComponent& C, StateComponent& R, StateComponent const& LeftHam, OperatorComponent const& HLeft, OperatorComponent const& HRight, StateComponent const& RightHam, PreExpansionAlgorithm Algo, int ExtraStates, int ExtraStatesPerSector, double RangeFindingOverhead)
{
   // quick return if we don't need to do anything
   if ((ExtraStates <= 0 && ExtraStatesPerSector <= 0) || Algo == PreExpansionAlgorithm::NoExpansion)
      return;

   if (Algo == PreExpansionAlgorithm::Random)
   {
      // Get the d*m dimensional full basis for the right hand side.  This is the candidate space for the states to add
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(R.LocalBasis(), R.Basis2()).Basis()).Basis();
      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(C.LocalBasis()), C.Basis1()).Basis()).Basis();

      // Get the dimensions of the 2-site null space
      // See comment for PreExpandBasis1
      // auto NumAvailablePerSector = subtract(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis2()));
      auto NumAvailablePerSector = subtract(CullMissingSectors(RFullBasis, LFullBasis), DimensionPerSector(C.Basis2()));

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis2()), ExtraStates, ExtraStatesPerSector, 1.0);

      //TRACE(ExpansionBasis.total_dimension())(RFullBasis.total_dimension()-C.Basis2().total_dimension())(R.Basis2().total_dimension());

      // Construct Haar random states to use as expansion vectors.
      MatrixOperator RMat = ReshapeBasis2(R);
      MatrixOperator M = MakeRandomMatrixOperator(ExpansionBasis, RMat.Basis2());
      MatrixOperator RMatNew = LQ_Factorize(tensor_col_sum(RMat, M)).second;

      // now we need to extend C to match
      MatrixOperator U = RMat * herm(RMatNew);
      C = C * U;

      R = ReshapeFromBasis2(RMatNew, R.LocalBasis(), R.Basis2());
   }
   else if (Algo == PreExpansionAlgorithm::RangeFinding || Algo == PreExpansionAlgorithm::RangeFindingProject)
   {
      // Get the d*m dimensional full basis for the right hand side.  This is the candidate space for the states to add
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(R.LocalBasis(), R.Basis2()).Basis()).Basis();
      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(C.LocalBasis()), C.Basis1()).Basis()).Basis();

      // Get the dimensions of the 2-site null space
      // See comment for PreExpandBasis1
      //auto NumAvailablePerSector = subtract(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis2()));
      auto NumAvailablePerSector = subtract(CullMissingSectors(RFullBasis, LFullBasis), DimensionPerSector(C.Basis2()));

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis2()), ExtraStates, ExtraStatesPerSector, 1.0);

      // Construct Gaussian random matrix
      StateComponent Omega = ReshapeFromBasis1(MakeRandomMatrixOperator(LFullBasis, ExpansionBasis), C.LocalBasis(), C.Basis1());

      // The RangeFindingProject algorithm projects out the kept states of C.
      if (Algo == PreExpansionAlgorithm::RangeFindingProject)
      {
         // In order to get the kept states of C, we can do a QR decomposition
         StateComponent LKeep = C;
         OrthogonalizeBasis2_QR(LKeep);
         Omega = Omega - LKeep * scalar_prod(herm(LKeep), Omega);
      }

      // Contract the left half
      StateComponent E = contract_from_left(HLeft, herm(Omega), LeftHam, C);
      StateComponent X = operator_prod_inner(HRight, E, R, herm(RightHam));
      // Project out the kept states from R
      X = X - scalar_prod(X, herm(R)) * R;
      // LQ decomposition, and we can throw away L
      OrthogonalizeBasis1_LQ(X);
      // Add the expansion vectors, and do another LQ factorization to ensure that the basis is orthogonal
      StateComponent RNew = tensor_col_sum(R, X);
      OrthogonalizeBasis1_LQ(RNew);

      // add the expansion vectors to C
      MatrixOperator U = scalar_prod(R, herm(RNew));
      C = C * U;
      R = RNew;
   }
   else
   {
      PANIC("Unsupported pre-expansion algorithm.");
   }
}

// Pre-expand Basis 1 of C, by adding vectors to the left-ortho Basis2() of L.
void
PreExpandBasis1(StateComponent& L, StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& HLeft, OperatorComponent const& HRight, StateComponent const& RightHam, PreExpansionAlgorithm Algo, int ExtraStates, int ExtraStatesPerSector, double RangeFindingOverhead)
{
   // quick return if we don't need to do anything
   if ((ExtraStates <= 0 && ExtraStatesPerSector <= 0) || Algo == PreExpansionAlgorithm::NoExpansion)
      return;

   if (Algo == PreExpansionAlgorithm::Random)
   {
      // Get the d*m dimensional full basis for the right hand side.  This is the candidate space for the states to add
      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(L.LocalBasis()), L.Basis1()).Basis()).Basis();
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(C.LocalBasis(), C.Basis2()).Basis()).Basis();

      // Get the dimensions of the 2-site null space, by subtracting the number of states in the kept basis from the candidate
      // space.  There is a choice here: we have LFullBasis states available,
      // but the maximum number of non-zero singular values is given by the min of the dimension of LFullBasis and RFullBasis.
      // Nevertheless we could keep more than the minimum, which will allow the truncation to mix between those states
      // when choosing the final kept states. But at minimum, we do want to cull states that have no available state at all
      // in RFullBasis.  Some limited empircal testing shows essentially no difference between these choices.
      // auto NumAvailablePerSector = subtract(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis1()));
      auto NumAvailablePerSector = subtract(CullMissingSectors(LFullBasis, RFullBasis), DimensionPerSector(C.Basis1()));

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis1()), ExtraStates, ExtraStatesPerSector, 1.0);

      //TRACE(ExpansionBasis.total_dimension())(RFullBasis.total_dimension()-C.Basis2().total_dimension())(R.Basis2().total_dimension());

      // Construct Haar random states to use as expansion vectors.
      MatrixOperator LMat = ReshapeBasis1(L);

      MatrixOperator M = MakeRandomMatrixOperator(LMat.Basis1(), ExpansionBasis);
      MatrixOperator LMatNew = QR_Factorize(tensor_row_sum(LMat, M)).first;

      // now we need to extend C to match
      MatrixOperator U = herm(LMatNew) * LMat;
      C = U * C;

      L = ReshapeFromBasis1(LMatNew, L.LocalBasis(), L.Basis1());
   }
   else if (Algo == PreExpansionAlgorithm::RangeFinding || Algo == PreExpansionAlgorithm::RangeFindingProject)
   {
      // Form the dm x k matrix, from doing a 2-site Hamiltonian matrix-vector multiply, with a dm x k matrix inserted

      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(L.LocalBasis()), L.Basis1()).Basis()).Basis();
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(C.LocalBasis(), C.Basis2()).Basis()).Basis();

      // See comment above
      //auto NumAvailablePerSector = subtract(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis1()));
      auto NumAvailablePerSector = subtract(CullMissingSectors(LFullBasis, RFullBasis), DimensionPerSector(C.Basis1()));

      // TODO: if we had rangefinding overhead...
      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis1()), ExtraStates, ExtraStatesPerSector, 1.0);

      // Construct Gaussian random matrix
      StateComponent Omega = ReshapeFromBasis2(MakeRandomMatrixOperator(ExpansionBasis, RFullBasis), C.LocalBasis(), C.Basis2());

      // The RangeFindingProject algorithm projects out the kept states of C.
      if (Algo == PreExpansionAlgorithm::RangeFindingProject)
      {
         // In order to get the kept states of C, we can do an LQ decomposition
         StateComponent RKeep = C;
         OrthogonalizeBasis1_LQ(RKeep);
         Omega = Omega - scalar_prod(Omega, herm(RKeep)) * RKeep;
      }

      // Contract the right half
      StateComponent F = contract_from_right(herm(HRight), Omega, RightHam, herm(C));
      StateComponent X = operator_prod_inner(HLeft, LeftHam, L, herm(F));
      // Project out the kept states from L
      X = X - L * scalar_prod(herm(L), X);
      // QR decomposition, we can throw away 'R'
      OrthogonalizeBasis2_QR(X);
      // Add the expansion vectors, and do another QR factorization to ensure that the basis is orthogonal
      StateComponent LNew = tensor_row_sum(L, X);
      OrthogonalizeBasis2_QR(LNew);

      // add the expansion vectors to C
      MatrixOperator U = scalar_prod(herm(LNew), L);
      C = U * C;
      L = LNew;
   }
   else
   {
      PANIC("Unsupported pre-expansion algorithm.");
   }
}


// Expand the Basis2 of C.
// The MPS is in the mixed canonical form cented around C.
// CRight must be right-orthonormal.
int
ExpandRightEnvironment(StateComponent& C, StateComponent& CRight,
                      StateComponent const& E, StateComponent const& F,
                      OperatorComponent const& HLeft, OperatorComponent const& HRight,
                      int StatesWanted, int ExtraStatesPerSector)
{
   #if 0
   int ExtraStates = StatesWanted - C.Basis2().total_dimension();

   // Calculate right null space of right site.
   StateComponent NRight = NullSpace1(CRight);

   // Perform SVD to left-orthogonalize the left site and extract singular value matrix.
   StateComponent COrtho = C;
   MatrixOperator VhLeft;
   RealDiagonalOperator DLeft;
   std::tie(DLeft, VhLeft) = OrthogonalizeBasis2(COrtho);

   // Now the 3S-inspired step: calculate X = new F matrix projected onto the null space.
   // Firstly project out the first and last columns of HRight.
   // NOTE: if the Hamiltonian is 2-dimensional MPO (i.e. there are no interactions here), then
   // the projector maps onto an empty set, so we need the 3-parameter constructor of the BasisList
   SimpleOperator Projector(BasisList(HRight.GetSymmetryList(), HRight.Basis1().begin()+1, HRight.Basis1().end()-1), HRight.Basis1());
   for (int i = 0; i < Projector.Basis1().size(); ++i)
      Projector(i, i+1) = 1.0;

   StateComponent X = contract_from_right(herm(Projector*HRight), NRight, F, herm(DLeft*VhLeft*CRight));
   StateComponent ELeft = contract_from_left(HLeft*herm(Projector), herm(COrtho), E, C);

   // Multiply each element of X by a prefactor depending on the corresponding element of E.
   for (int i = 0; i < X.size(); ++i)
   {
      // We add an epsilon term to the prefactor to avoid the corner case where
      // the prefactor is zero for all elements, and any possible new states
      // are ignored.
      double Prefactor = norm_frob(ELeft[i]) + PrefactorEpsilon;
      X[i] *= Prefactor;
   }

   // Take the truncated SVD of X.  This appears asymmetric compared with ExpandLeftEnvironment, but it
   // is actually OK: the equivalent 'reflected' operation would be to ReshapeBasis1(herm(X)), but instead
   // we can just ReshapeBasis2(X), so the rest of the code is essentially identical to ExpandLeftEnvironment,
   // except for swapping C/CRight and Basis1/Basis2
   MatrixOperator XExpand = ReshapeBasis2(X);
   XExpand = MatrixOperator(XExpand.Basis1(), VectorBasis(XExpand.GetSymmetryList()));

   CMatSVD SVD(XExpand, CMatSVD::Left);

   auto StatesToKeep = TruncateExtraStates(SVD.begin(), SVD.end(), ExtraStates, ExtraStatesPerSector, true);

   MatrixOperator U = SVD.ConstructLeftVectors(StatesToKeep.begin(), StatesToKeep.end());

   // Construct new basis.
   SumBasis<VectorBasis> NewBasis(CRight.Basis1(), U.Basis2());
   // Regularize the new basis.
   Regularizer R(NewBasis);

   // Add the new states to CRight, and add zeros to C.
   CRight = RegularizeBasis1(R, tensor_col_sum(CRight, prod(herm(U), NRight), NewBasis));

   StateComponent Z = StateComponent(C.LocalBasis(), C.Basis1(), U.Basis2());
   C = RegularizeBasis2(tensor_row_sum(C, Z, NewBasis), R);

   return StatesToKeep.size();
   #endif
   return 0;
}

#if 0
void DMRG::ExpandEnvironmentLeft(int ExtraStates)
{
   // Expand the environment for Basis2 of *C

   // Merge E and A matrices
   StateComponent A = *C;
   StateComponent EA = local_tensor_prod(LeftHamiltonian.back(), A);
   // EA is m x (w*d) x m

   ProductBasis<BasisList, BasisList> PBasis(A.LocalBasis(), EA.LocalBasis());

   SimpleOperator Ham = reshape_to_matrix(*H);

   // Apply the MPO
   EA = local_prod(herm(Ham), EA);
   // EA is now m x (d*w) x m

   // project onto null space
   EA -= local_tensor_prod(A, contract_local_tensor_prod_left(herm(A), EA, PBasis));

   // Pre-SVD to make the final SVD faster
   TruncateBasis2(EA);

   // SVD again to get the expanded basis
   AMatSVD SL(EA, PBasis);

   // How to choose which states to keep?
   // If we do any expansion at all then we ought to keep at least one state in each possible quantum number sector, if possible.
   TruncationInfo Info;
   RealDiagonalOperator L;
   StateComponent X;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfo, Info);
   SL.ConstructMatrices(SL.begin(), Cutoff, A, L, X);  // L and X are not used here

   // A is now the expanded set of states.  Merge it into *C
   A = tensor_row_sum(*C, A, SumBasis<VectorBasis>(C->Basis2(), A.Basis2()));

   // The added rows are not going to be exactly orthogonal to the existing states, especially if we forced a few states that
   // have zero singular values. We ought to do another SVD, or a QR decomposition to orthogonalize the basis.

   // Reconstruct the Hamiltonian
   // In principle we could just construct the new rows using A, but easier to just reconstruct the full matrices.

   // We also need Lambda in the new basis


}
#endif

int DMRG::ExpandLeftEnvironment(int StatesWanted, int ExtraStatesPerSector)
{
   auto CLeft = C;
   --CLeft;
   HamMatrices.pop_left();
   auto HLeft = H;
   --HLeft;

   //::ExpandLeftEnvironment(*CLeft, *C, HamMatrices.left(), HamMatrices.right(), *HLeft, *H, StatesWanted, ExtraStatesPerSector);

   PreExpandBasis1(*CLeft, *C, HamMatrices.left(), *HLeft, *H, HamMatrices.right(), PreExpansionAlgo, StatesWanted-C->Basis1().total_dimension(), ExtraStatesPerSector, 1.0);

   // reconstruct the E matrix with the expanded environment
   HamMatrices.push_left(contract_from_left(*HLeft, herm(*CLeft), HamMatrices.left(), *CLeft));

   return C->Basis1().total_dimension();
}

int DMRG::ExpandRightEnvironment(int StatesWanted, int ExtraStatesPerSector)
{
   auto CRight = C;
   ++CRight;
   HamMatrices.pop_right();
   auto HRight = H;
   ++HRight;

   PreExpandBasis2(*C, *CRight, HamMatrices.left(), *H, *HRight, HamMatrices.right(), PreExpansionAlgo, StatesWanted-C->Basis2().total_dimension(), ExtraStatesPerSector, 1.0);

   // auto EnvStates = ::ExpandRightEnvironment(*C, *CRight, HamMatrices.left(), HamMatrices.right(), *H, *HRight, StatesWanted, ExtraStatesPerSector);

   // reconstsruct the F matrix with the expanded environment
   HamMatrices.push_right(contract_from_right(herm(*HRight), *CRight, HamMatrices.right(), herm(*CRight)));

   return C->Basis2().total_dimension();
}

TruncationInfo DMRG::TruncateAndShiftLeft(StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector)
{
   TruncationInfo Info;
   MatrixOperator X = TruncateExpandBasis1(*C, HamMatrices.left(), *H, HamMatrices.right(), PostExpansionAlgo, States, ExtraStates, ExtraStatesPerSector, Info, RangeFindingOverhead);
   if (Verbose > 1)
   {
      std::cerr << "Truncating left basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   HamMatrices.push_right(contract_from_right(herm(*H), *C, HamMatrices.right(), herm(*C)));
   HamMatrices.pop_left();

   // next site
   --Site;
   --H;
   --C;

   *C = prod(*C, X);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

TruncationInfo DMRG::TruncateAndShiftRight(StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector)
{
   // Truncate right
   StateComponent CC = *C;
   TruncationInfo Info;
   MatrixOperator X = TruncateExpandBasis2(CC, HamMatrices.left(), *H, HamMatrices.right(), PostExpansionAlgo, States, ExtraStates, ExtraStatesPerSector, Info, RangeFindingOverhead);
   *C = CC;
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }
   // update blocks
   HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));
   HamMatrices.pop_right();

   // next site
   ++Site;
   ++H;
   ++C;

   *C = prod(X, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

void DMRG::debug_check_structure() const
{
#if 0
   DEBUG_CHECK_EQUAL(HamMatrices.Left().Basis1(), Psi.Center().Basis1());
   DEBUG_CHECK_EQUAL(HamMatrices.Right().Basis1(), Psi.Center().Basis2());

   DEBUG_CHECK_EQUAL(HamMatrices.Left().Basis2(), Psi.Center().Basis1());
   DEBUG_CHECK_EQUAL(HamMatrices.Right().Basis2(), Psi.Center().Basis2());

   for (std::size_t j = 0; j < Ortho.size(); ++j)
   {
      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Left().Basis1(), Psi.Center().Basis1());
      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Right().Basis1(), Psi.Center().Basis2());

      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Left().Basis2(), Ortho[j].Center().Basis1());
      DEBUG_CHECK_EQUAL(PsiOrthoProjector[j].Right().Basis2(), Ortho[j].Center().Basis2());
   }
#endif
}
