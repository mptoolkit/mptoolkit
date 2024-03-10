// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/dmrg.cpp
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
constexpr std::array<char const*,5> PreExpansionTraits::Names;
constexpr std::array<char const*, 6> PostExpansionTraits::Names;

std::ostream& operator<<(std::ostream& out, OversamplingInfo const& k)
{
   out << "Oversampling min: " << k.Add << " factor: " << k.Scale << " per-sector: " << k.ExtraPerSector << '\n';
   return out;
}

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

// helper to exclude y from x.  That is, the elements of Result' satisfy
// Result'[a] = std::max(x[a] - y[a], 0)
// This is used in calculating the dimensions of sectors, where we want to exclude (for example) states that are already kept
template <typename T>
std::map<T, int>
exclude(std::map<T, int> x, std::map<T, int> const& y)
{
   for (auto const& yi : y)
   {
      if ((x[yi.first] -= yi.second) <= 0)
         x.erase(yi.first);
   }
   return x;
}

// Get a basis for use in the randomized SVD, selecting ExtraStates*RangeFindingOverhead states out of the available states.
// The quantum number sectors for the singular vectors are chosen according to the relative number of states in
// KeptStateDimension, with at least ExtraStatesPerSector states in each available quantum number sector.
// NumAvailablePerSector is the total size of the space; we exclude the KeptStateDimension before allocating states.
// That is, the number of states allocated to some quantum number sector q is no larger than
// NumAvailablePerSector[q] - KeptStateDimension[q]
// Typically KeptStateDimension will be initialized as DimensionPerSector() of some basis of kept states.
// On entry: X is a p x q matrix (typically = p m*d, q = w*m)
// On exit: Result' is a VectorBasis of size ExtraStates + additional states per sector
VectorBasis
MakeExpansionBasis(QuantumNumbers::SymmetryList const& SL, std::map<QuantumNumbers::QuantumNumber, int> NumAvailablePerSector, std::map<QuantumNumber, int> const& KeptStateDimension, int ExtraStates, int ExtraStatesPerSector, OversamplingInfo const& Oversampling)
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
   NumAvailablePerSector = exclude(NumAvailablePerSector, KeptStateDimension);
   auto NumExtraPerSector = DistributeStates(Weights, NumAvailablePerSector, ExtraStates, ExtraStatesPerSector+Oversampling.ExtraPerSector);

   // Now that we have the distribution, apply the over-sampling
   for (auto& r : NumExtraPerSector)
   {
      r.second = std::min(NumAvailablePerSector.at(r.first), Oversampling(r.second));
   }
   return VectorBasis(SL, NumExtraPerSector.begin(), NumExtraPerSector.end());
}

#if 0
// The randomized SVD version of the 3S algorithm
MatrixOperator
TruncateMixBasis1(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int k, OversamplingInfo const& Oversampling, TruncationInfo& Info)
{
   // Reshape C to m x dm
   MatrixOperator CMat = ReshapeBasis2(C);

   CMatSVD DM(CMat);
   // unless we are reducing the basis size, this will keep all of the non-sero singular values.  But we need to
   // do the SVD anyway so that we have VKeep, of the states to keep.
   auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
   MatrixOperator VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
   StateComponent R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());
   // R is the right-orthogonal A-matrix of states to keep.  The algorithms below obtain R in different ways.

   SimpleOperator P = ProjectRemoveFirstLast(LeftHam.LocalBasis());
   SimpleOperator W = ConstructWeights(local_prod(P, LeftHam));

   // Construct the embedding matrix from the w*m dimensional basis to k dimensions. Firstly get the w*m basis
   VectorBasis RBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(W.Basis1(), C.Basis1()).Basis()).Basis();

   // Get the basis of candidate states based on the ranks of the kept states
   VectorBasis EBasis = MakeExpansionBasis(VKeep.GetSymmetryList(), RankPerSector(VKeep.Basis2(), RBasis), DimensionPerSector(VKeep.Basis1()), k, Oversampling);

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
   // We don't have to worry about RNull containing numerically small elements in the direction of R, since
   // we're going to combine it with the C tensor anyway

   // Now find the matrix elements of the mixing term in the RNull basis
   StateComponent D = contract_from_right(herm(W*P*H), RNull, RightHam, herm(C));
   MatrixOperator X = ReshapeBasis2(D);

   // Now mix X into C and do another SVD of the combined set

   // normalize X
   // We choose the mixing factor from the truncation error associated with the final k singular values of C.
   // If the final k singular values are all zero, then we are going to keep all of the expansion vectors anyway so
   // the choice of normalization is not sensitive; we can set it to some small number such as 1e-16 (since it
   // is used in an SVD we take the square root, so numerically epsilon ends ip as 1e-8)
   // Pivot marks the start of the last k singular values.  But clamp it between 0 and the number of kept states
   int Pivot = std::min(std::max(C.Basis1().total_dimension() - k, 0), R.Basis1().total_dimension());
   double MixFactor = std::max(DensitySigma(DM.begin()+Pivot, DM.end(), 1e-16);
   X *= (std::sqrt(MixFactor) / norm_frob(X));

   // Now reconstruct R using the SVD of the mixing matrix
   DM = CMatSVD(tensor_row_sum(CMat, X));
   DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
   VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
   R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());

   // Construct the new Lambda matrix
   MatrixOperator Lambda = scalar_prod(C, herm(R));
   C = std::move(R);
   CHECK_EQUAL(LeftHam.Basis1(), Lambda.Basis1());
   CHECK_EQUAL(Lambda.Basis2(), C.Basis1());
   CHECK_EQUAL(C.Basis2(), RightHam.Basis1());
   return Lambda;
}
#endif

// Expand the Basis1 of C.
// On exit, Result' * C' = C (up to truncation!), and C is right-orthogonal
MatrixOperator
TruncateExpandBasis1(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector, TruncationInfo& Info, OversamplingInfo Oversampling)
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

   if (Algo == PostExpansionAlgorithm::NoExpansion)
   {
      CMatSVD DM(CMat);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
      R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());
      Info.ExtraStates_ = 0;
   }
   else if (Algo == PostExpansionAlgorithm::RSVD || Algo == PostExpansionAlgorithm::RangeFinding)
   {
      // We can do a thin SVD here.  In the no expansion case we just do the basic truncation and we are done.
      // In the fastrangefinding algorithm we implicitly use the full space of m*d states, but project out
      // the kept states.  The arrangement of contractions means that this is faster than constructing the
      // null space explicitly.
      CMatSVD DM(CMat);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
      R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());

      SimpleOperator P = ProjectRemoveFirstLast(LeftHam.LocalBasis());
      SimpleOperator W = ConstructWeights(local_prod(P, LeftHam));

      // Construct the embedding matrix from the w*m dimensional basis to k dimensions. Firstly get the w*m basis
      VectorBasis RBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(W.Basis1(), C.Basis1()).Basis()).Basis();

      // Get the basis of candidate states based on the ranks of the kept states
      if (Algo == PostExpansionAlgorithm::RangeFinding)
         Oversampling = OversamplingInfo();  // Don't oversample in the rangefinding algorithm

      VectorBasis EBasis = MakeExpansionBasis(VKeep.GetSymmetryList(), RankPerSector(VKeep.Basis2(), RBasis), DimensionPerSector(VKeep.Basis1()), ExtraStates, ExtraStatesPerSector, Oversampling);

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

      if (Algo == PostExpansionAlgorithm::RSVD)
      {
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
      }
      Info.ExtraStates_ = RNull.Basis1().total_dimension();

      R = RegularizeBasis1(tensor_col_sum(R, RNull));
      OrthogonalizeBasis1_LQ(R);  // additional orthogonalization step, ensure vectors are all orthogonal
   }
   else if (Algo == PostExpansionAlgorithm::SVD || Algo == PostExpansionAlgorithm::Random)
   {
      // The SVD, rangefinding, and random algorithms differ only by rangefinding using a randomized SVD.
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

         CMatSVD ExpandDM(X, CMatSVD::Left);
         auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
         QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
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
TruncateExpandBasis2(StateComponent& C, StateComponent const& LeftHam, OperatorComponent const& H, StateComponent const& RightHam, PostExpansionAlgorithm Algo, StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector, TruncationInfo& Info, OversamplingInfo Oversampling)
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

   if (Algo == PostExpansionAlgorithm::NoExpansion)
   {
      CMatSVD DM(CMat);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator UKeep = DM.ConstructLeftVectors(DM.begin(), DMPivot);
      L = ReshapeFromBasis1(UKeep, C.LocalBasis(), C.Basis1());
      Info.ExtraStates_ = 0;
   }
   else if (Algo == PostExpansionAlgorithm::RSVD || Algo == PostExpansionAlgorithm::RangeFinding)
   {
      // if we're range-finding, reset the Oversampling to zero
      if (Algo == PostExpansionAlgorithm::RangeFinding)
         Oversampling = OversamplingInfo();
      // We can do a thin SVD here.  In the no expansion case we just do the basic truncation and we are done.
      // In the fastrangefinding algorithm we implicitly use the full space of m*d states, but project out
      // the kept states.  The arrangement of contractions means that this is faster than constructing the
      // null space explicitly.
      CMatSVD DM(CMat);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator UKeep = DM.ConstructLeftVectors(DM.begin(), DMPivot);
      L = ReshapeFromBasis1(UKeep, C.LocalBasis(), C.Basis1());

      SimpleOperator P = ProjectRemoveFirstLast(RightHam.LocalBasis());
      SimpleOperator W = ConstructWeights(local_prod(P, RightHam));

      // Construct the embedding matrix from the w*m dimensional basis to k dimensions. Firstly get the w*m basis
      VectorBasis RBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(W.Basis2(), C.Basis2()).Basis()).Basis();

      // Get the basis of candidate states based on the ranks of the kept states
      if (Algo == PostExpansionAlgorithm::RangeFinding)
         Oversampling = OversamplingInfo();  // Don't oversample in the rangefinding algorithm

      VectorBasis EBasis = MakeExpansionBasis(UKeep.GetSymmetryList(), RankPerSector(UKeep.Basis1(), RBasis), DimensionPerSector(UKeep.Basis2()), ExtraStates, ExtraStatesPerSector, Oversampling);

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

      // At this point, LNull might have some non-zero numerical overlap with L. But it should be very small, since
      // the number of expansion vectors is controlled by the rank of the null space.
      // Since we do a QR at the end anyway, we don't need to worry about a numerically non-orthogonal basis.
      // A faster alternative to the QR would be Gram-Schmidt orthogonalization.

      if (Algo == PostExpansionAlgorithm::RSVD)
      {
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
      }

      Info.ExtraStates_ = LNull.Basis2().total_dimension();

      // LNull might not be exactly orthogonal
      L = RegularizeBasis2(tensor_row_sum(L, LNull));
      OrthogonalizeBasis2_QR(L);
   }
   else if (Algo == PostExpansionAlgorithm::SVD || Algo == PostExpansionAlgorithm::Random)
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
         // For SVD, we construct the mixing term from the 3S algorithm
         SimpleOperator P = ProjectRemoveFirstLast(RightHam.LocalBasis());
         SimpleOperator W = ConstructWeights(local_prod(P, RightHam));
         MatrixOperator X = ReshapeBasis2(contract_from_left(H*herm(P)*W, herm(LNull), LeftHam, C));
         // X is now (discarded basis) x (wm) matrix

         // Get the kept and discarded states
         CMatSVD ExpandDM(X, CMatSVD::Left);
         auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
         QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
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
     Oversampling(),
     ProjectTwoSiteTangent(false)
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
   return std::numeric_limits<double>::signaling_NaN();
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

   // this is where we would write some log info per sweep
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

std::pair<std::complex<double>, TruncationInfo>
DMRG::SolveCoarseGrainRight(StatesInfo const& States)
{
   auto R = C; ++R;
   auto CoarseGrainBasis = make_product_basis(C->LocalBasis(), R->LocalBasis());
   StateComponent X = local_tensor_prod(*C, *R);

   auto HR = H; ++HR;
   OperatorComponent HX = local_tensor_prod(*H, *HR);

   HamMatrices.pop_right();

   Solver_.Solve(X, HamMatrices.left(), HX, HamMatrices.right());

   // NOTE: *C is typically nearly, but not exactly normalized at this point.
   // There is no point normalizing here, since we need to do it anyway after the truncation.

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   IterationEnergyBeforeTrunc = IterationEnergy;

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   IterationEnergyBeforeTrunc = IterationEnergy;

   AMatSVD DM(X, CoarseGrainBasis);

   TruncationInfo Info;
   auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
   RealDiagonalOperator Lambda;
   DM.ConstructMatrices(DM.begin(), DMPivot, *C, Lambda, *R);

   Lambda *= 1.0 / norm_frob(Lambda); // normalize

   // update block Hamiltonian
   HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));

   // next site
   ++Site;      // note: this is now one site ahead of the site that was optimized.
   ++C;
   *C = Lambda * (*C);
   ++H;

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return std::make_pair(IterationEnergy, Info);
}

std::pair<std::complex<double>, TruncationInfo>
DMRG::SolveCoarseGrainLeft(StatesInfo const& States)
{

   auto L = C; --L;
   auto CoarseGrainBasis = make_product_basis(L->LocalBasis(), C->LocalBasis());
   StateComponent X = local_tensor_prod(*L, *C);

   auto HL = H; --HL;
   OperatorComponent HX = local_tensor_prod(*HL, *H);

   HamMatrices.pop_left();

   Solver_.Solve(X, HamMatrices.left(), HX, HamMatrices.right());

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   IterationEnergyBeforeTrunc = IterationEnergy;

   AMatSVD DM(X, CoarseGrainBasis);

   TruncationInfo Info;
   auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
   RealDiagonalOperator Lambda;
   DM.ConstructMatrices(DM.begin(), DMPivot, *L, Lambda, *C);

   Lambda *= 1.0 / norm_frob(Lambda); // normalize

   // update block Hamiltonian
   HamMatrices.push_right(contract_from_right(herm(*H), *C, HamMatrices.right(), herm(*C)));

   // next site
   --Site;         // note: this is now one site ahead of the site that was optimized.
   --C;
   *C = (*C) * Lambda;
   --H;
   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return std::make_pair(IterationEnergy, Info);
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

#if 0
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
}
#endif

void OrthogonalizeRowsAgainst(StateComponent& X, StateComponent const& Y)
{
   #if 1
   // This is faster than reshaping
   X = X - scalar_prod(X, herm(Y)) * Y;
   OrthogonalizeBasis1_LQ(X);
   X = X - scalar_prod(X, herm(Y)) * Y;
   OrthogonalizeBasis1_LQ(X);
   #else
   auto XM = ReshapeBasis2(X);
   auto YM = ReshapeBasis2(Y);
   OrthogonalizeRowsAgainst(XM, YM);
   X = ReshapeFromBasis2(XM, X.LocalBasis(), X.Basis2());
   #endif
}

void OrthogonalizeColsAgainst(StateComponent& X, StateComponent const& Y)
{
   #if 1
   // This is faster than reshaping
   X = X - Y * scalar_prod(herm(Y), X);
   OrthogonalizeBasis2_QR(X);
   X = X - Y * scalar_prod(herm(Y), X);
   OrthogonalizeBasis2_QR(X);
   #else
   auto XM = ReshapeBasis1(X);
   auto YM = ReshapeBasis1(Y);
   OrthogonalizeColsAgainst(XM, YM);
   X = ReshapeFromBasis1(XM, X.LocalBasis(), X.Basis1());
   #endif
}

// Pre-expand Basis 2 of C, by adding vectors to the right-ortho Basis1() of R.
// Returns the expansion vectors in the basis of R.
// Result'.Basis2() == R.Basis2()
// Result'.Basis1() can be anything that doesn't exhaust the rank of R.Basis1() + Result'.Basis1(),
// eg returning random vectors is OK, they do not need to be in the null space as it is assumed that the
// caller will handle that.
StateComponent
PreExpandBasis2(StateComponent const& C, StateComponent const& R, StateComponent const& LeftHam, OperatorComponent const& HLeft, OperatorComponent const& HRight, StateComponent const& RightHam, PreExpansionAlgorithm Algo, int ExtraStates, int ExtraStatesPerSector, OversamplingInfo Oversampling, bool ProjectTwoSiteTangent)
{
   // quick return if we don't need to do anything
   if ((ExtraStates <= 0 && ExtraStatesPerSector <= 0) || Algo == PreExpansionAlgorithm::NoExpansion)
      return StateComponent(R.LocalBasis(), VectorBasis(R.GetSymmetryList()), R.Basis2());

   if (Algo == PreExpansionAlgorithm::SVD)
   {
      StateComponent RNull = NullSpace1(R);
      StateComponent F = contract_from_right(herm(HRight), RNull, RightHam, herm(R));
      StateComponent CPrime = operator_prod_inner(HLeft, LeftHam, C, herm(F));

      if (ProjectTwoSiteTangent)
      {
         StateComponent LKeep = C;
         OrthogonalizeBasis2_QR(LKeep);
         CPrime = CPrime - LKeep * scalar_prod(herm(LKeep), CPrime);
      }
      MatrixOperator X = ReshapeBasis1(CPrime);

      CMatSVD ExpandDM(X, CMatSVD::Right);
      auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, true);
      auto QExpand = ExpandDM.ConstructRightVectors(ExpandedStates.begin(), ExpandedStates.end());

      StateComponent Q = QExpand * RNull;
      return Q;
   }

   else if (Algo == PreExpansionAlgorithm::Random)
   {
      // Get the d*m dimensional full basis for the right hand side.  This is the candidate space for the states to add
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(R.LocalBasis(), R.Basis2()).Basis()).Basis();
      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(C.LocalBasis()), C.Basis1()).Basis()).Basis();

      // Get the dimensions of the 2-site null space
      // See comment for PreExpandBasis1
      // auto NumAvailablePerSector = exclude(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis2()));
      auto NumAvailablePerSector = CullMissingSectors(RFullBasis, LFullBasis);

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis2()), ExtraStates, ExtraStatesPerSector, OversamplingInfo());  // don't oversample

      // Construct Haar random states to use as expansion vectors.
      MatrixOperator RMat = ReshapeBasis2(R);
      MatrixOperator M = MakeRandomMatrixOperator(ExpansionBasis, RMat.Basis2());

      return ReshapeFromBasis2(M, R.LocalBasis(), R.Basis2());
   }
   else if (Algo == PreExpansionAlgorithm::RangeFinding || Algo == PreExpansionAlgorithm::RSVD)
   {
      if (Algo == PreExpansionAlgorithm::RangeFinding)
         Oversampling = OversamplingInfo();

      // Get the d*m dimensional full basis for the right hand side.  This is the candidate space for the states to add
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(R.LocalBasis(), R.Basis2()).Basis()).Basis();
      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(C.LocalBasis()), C.Basis1()).Basis()).Basis();

      // Get the dimensions of the 2-site null space
      // See comment for PreExpandBasis1
      //auto NumAvailablePerSector = exclude(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis2()));
      auto NumAvailablePerSector = CullMissingSectors(RFullBasis, LFullBasis);

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis2()), ExtraStates, ExtraStatesPerSector, Oversampling);

      // Construct Gaussian random matrix
      StateComponent Omega = ReshapeFromBasis1(MakeRandomMatrixOperator(LFullBasis, ExpansionBasis), C.LocalBasis(), C.Basis1());

      // The RangeFindingProject algorithm projects out the kept states of C.
      if (ProjectTwoSiteTangent)
      {
         // In order to get the kept states of C, we need to orthogonalize it, and we can do a QR decomposition
         StateComponent LKeep = C;
         OrthogonalizeBasis2_QR(LKeep);
         Omega = Omega - LKeep * scalar_prod(herm(LKeep), Omega);
      }

      // Contract the left half
      StateComponent E = contract_from_left(HLeft, herm(Omega), LeftHam, C);
      StateComponent Q = operator_prod_inner(HRight, E, R, herm(RightHam));

      if (Algo == PreExpansionAlgorithm::RSVD)
      {
         // Complete the range finding step and feed Q^\dagger back into the tensor network
         //OrthogonalizeRowsAgainst(Q, R)
         Q = Q - scalar_prod(Q, herm(R)) * R;
         OrthogonalizeBasis1_LQ(Q);

         StateComponent F = contract_from_right(herm(HRight), Q, RightHam, herm(R));
         StateComponent QX = operator_prod_inner(HLeft, LeftHam, C, herm(F));

         MatrixOperator QXf = ReshapeBasis1(QX);
         CMatSVD ExpandDM(QXf, CMatSVD::Right);
         auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, true);
         auto QExpand = ExpandDM.ConstructRightVectors(ExpandedStates.begin(), ExpandedStates.end());
         Q = QExpand*Q;
      }

      return Q;
   }
   else
   {
      PANIC("Unsupported pre-expansion algorithm.");
      return C; // suppress warning about missing return
   }
}

// Pre-expand Basis 1 of C, by adding vectors to the left-ortho Basis2() of L.
// Returns the expansion vectors in the basis of L.
// Result'.Basis1() == L.Basis1()
// Result'.Basis2() can be anything that doesn't exhaust the rank of L.Basis2() + Result'.Basis2(),
// eg returning random vectors is OK, they do not need to be in the null space as it is assumed that the
// caller will handle that.
StateComponent
PreExpandBasis1(StateComponent const& L, StateComponent const& C, StateComponent const& LeftHam, OperatorComponent const& HLeft, OperatorComponent const& HRight, StateComponent const& RightHam, PreExpansionAlgorithm Algo, int ExtraStates, int ExtraStatesPerSector, OversamplingInfo Oversampling, bool ProjectTwoSiteTangent)
{
   // quick return if we don't need to do anything
   if ((ExtraStates <= 0 && ExtraStatesPerSector <= 0) || Algo == PreExpansionAlgorithm::NoExpansion)
      return StateComponent(L.LocalBasis(), L.Basis1(), VectorBasis(L.GetSymmetryList()));

   if (Algo == PreExpansionAlgorithm::SVD)
   {
      // This is slow, we could avoid constructing the null space but the bottleneck is the SVD anyway
      StateComponent LNull = NullSpace2(L);
      StateComponent E = contract_from_left(HLeft, herm(LNull), LeftHam, L);
      StateComponent CPrime = operator_prod_inner(HRight, E, C, herm(RightHam));

      if (ProjectTwoSiteTangent)
      {
         StateComponent RKeep = C;
         OrthogonalizeBasis1_LQ(RKeep);
         CPrime = CPrime - scalar_prod(CPrime, herm(RKeep)) * RKeep;
      }
      MatrixOperator X = ReshapeBasis2(CPrime);

      CMatSVD ExpandDM(X, CMatSVD::Left);
      auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, true);
      auto QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());

      StateComponent Q = LNull * QExpand;
      return Q;
   }
   else if (Algo == PreExpansionAlgorithm::Random)
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
      // auto NumAvailablePerSector = exclude(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis1()));
      auto NumAvailablePerSector = CullMissingSectors(LFullBasis, RFullBasis);

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis1()), ExtraStates, ExtraStatesPerSector, OversamplingInfo()); // Don't oversample

      // Construct Haar random states to use as expansion vectors.
      MatrixOperator LMat = ReshapeBasis1(L);
      MatrixOperator M = MakeRandomMatrixOperator(LMat.Basis1(), ExpansionBasis);

      return ReshapeFromBasis1(M, L.LocalBasis(), L.Basis1());
   }
   else if (Algo == PreExpansionAlgorithm::RangeFinding || Algo == PreExpansionAlgorithm::RSVD)
   {
      if (Algo == PreExpansionAlgorithm::RangeFinding)
         Oversampling = OversamplingInfo();

      // Form the dm x k matrix, from doing a 2-site Hamiltonian matrix-vector multiply, with a dm x k matrix inserted

      VectorBasis LFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(adjoint(L.LocalBasis()), L.Basis1()).Basis()).Basis();
      VectorBasis RFullBasis = Regularizer(ProductBasis<BasisList, VectorBasis>(C.LocalBasis(), C.Basis2()).Basis()).Basis();

      // See comment above
      //auto NumAvailablePerSector = exclude(RankPerSector(LFullBasis, RFullBasis), DimensionPerSector(C.Basis1()));
      // Note that we no longer need to exclude DimensionPerSector(C.Basis1()) from NumAvailablePerSector, since this
      // is now done by MakeExpansionBasis()
      auto NumAvailablePerSector = CullMissingSectors(LFullBasis, RFullBasis);

      VectorBasis ExpansionBasis = MakeExpansionBasis(C.GetSymmetryList(), NumAvailablePerSector, DimensionPerSector(C.Basis1()), ExtraStates, ExtraStatesPerSector, Oversampling);

      // Construct Gaussian random matrix
      StateComponent Omega = ReshapeFromBasis2(MakeRandomMatrixOperator(ExpansionBasis, RFullBasis), C.LocalBasis(), C.Basis2());

      if (ProjectTwoSiteTangent)
      {
         // In order to get the kept states of C, we need to orthogonalize it, and we can do an LQ decomposition
         StateComponent RKeep = C;
         OrthogonalizeBasis1_LQ(RKeep);
         Omega = Omega - scalar_prod(Omega, herm(RKeep)) * RKeep;
      }

      // Contract the right half
      StateComponent F = contract_from_right(herm(HRight), Omega, RightHam, herm(C));
      StateComponent Q = operator_prod_inner(HLeft, LeftHam, L, herm(F));

      if (Algo == PreExpansionAlgorithm::RSVD)
      {
         // Complete the range finding step and feed Q^\dagger back into the tensor network
         //OrthogonalizeColsAgainst(Q, L);
         Q = Q - L * scalar_prod(herm(L), Q);
         OrthogonalizeBasis2_QR(Q);

         StateComponent E = contract_from_left(HLeft, herm(Q), LeftHam, L);
         StateComponent QX = operator_prod_inner(HRight, E, C, herm(RightHam));

         MatrixOperator QXf = ReshapeBasis2(QX);
         CMatSVD ExpandDM(QXf, CMatSVD::Left);
         auto ExpandedStates = TruncateExtraStates(ExpandDM.begin(), ExpandDM.end(), ExtraStates, ExtraStatesPerSector, true);
         auto QExpand = ExpandDM.ConstructLeftVectors(ExpandedStates.begin(), ExpandedStates.end());
         Q = Q*QExpand;
      }

      return Q;
   }
   else
   {
      PANIC("Unsupported pre-expansion algorithm.");
      return C; // suppress warning about missing return
   }
}

int DMRG::ExpandLeftEnvironment(int StatesWanted, int ExtraStatesPerSector)
{
   auto CLeft = C;
   --CLeft;
   StateComponent E = HamMatrices.left(); // the old E matrices
   HamMatrices.pop_left();
   auto HLeft = H;
   --HLeft;

   StateComponent CC = *C;
   StateComponent L = *CLeft;

   StateComponent LExpand = PreExpandBasis1(L, CC, HamMatrices.left(), *HLeft, *H, HamMatrices.right(), PreExpansionAlgo, StatesWanted-C->Basis1().total_dimension(), ExtraStatesPerSector, Oversampling, ProjectTwoSiteTangent);

   #if 0
   // Orthogonalize the expansion vectors and add them to the left side
   OrthogonalizeColsAgainst(LExpand, L);
   StateComponent LNew = tensor_row_sum(L, LExpand);
   *CLeft = LNew;

   // Augment the C matrix with zeros
   StateComponent CExpand(CC.LocalBasis(), LExpand.Basis2(), CC.Basis2());
   *C = tensor_col_sum(CC, CExpand);

   // Expansion of the E matrices
   // There are 3 components we need to add.  The additional rows for (LExpand, CLeft),
   // the additional columns for (CLeft, LExpand), and the diagonal (LExpand, LExpand).
   // We can add the rows first, and then add the columns (including the diagonal)
   // Alternatively, we can reconstruct the new E matrices from LNew
   StateComponent ERows = contract_from_left(*HLeft, herm(LExpand), HamMatrices.left(), L);
   E = tensor_col_sum(E, ERows);
   StateComponent ECols = contract_from_left(*HLeft, herm(LNew), HamMatrices.left(), LExpand);
   E = tensor_row_sum(E, ECols);
   HamMatrices.push_left(E);

   #else
   StateComponent LNew = tensor_row_sum(L, LExpand);
   OrthogonalizeBasis2_QR(LNew);

   // The full reconstruction would be:
   HamMatrices.push_left(contract_from_left(*HLeft, herm(LNew), HamMatrices.left(), LNew));

   // Update C in the new basis
   MatrixOperator U = scalar_prod(herm(LNew), L);
   *C = U * CC;
   *CLeft = LNew;

   #endif

   return C->Basis1().total_dimension();
}

int DMRG::ExpandRightEnvironment(int StatesWanted, int ExtraStatesPerSector)
{
   auto CRight = C;
   ++CRight;
   StateComponent F = HamMatrices.right();  // the old F matrix
   HamMatrices.pop_right();
   auto HRight = H;
   ++HRight;

   StateComponent CC = *C;
   StateComponent R = *CRight;

   StateComponent RExpand = PreExpandBasis2(CC, R, HamMatrices.left(), *H, *HRight, HamMatrices.right(), PreExpansionAlgo, StatesWanted-C->Basis2().total_dimension(), ExtraStatesPerSector, Oversampling, ProjectTwoSiteTangent);

   #if 0
   // In theory this code should be faster, but somehow it isn't. Enabling this makes the code slower,
   // and the CPU time per step is rather sporadic.

   // Orthogonalize the expansion vectors and add add them as new rows
   OrthogonalizeRowsAgainst(RExpand, R);
   StateComponent RNew = tensor_col_sum(R, RExpand);
   *CRight = RNew;

   // Augment the C matrix with zeros
   StateComponent CExpand(CC.LocalBasis(), CC.Basis1(), RExpand.Basis1());
   *C = tensor_row_sum(CC, CExpand);

   // Expansion of the F matrices
   // There are 3 components we need to add.  The additional rows for (RExpand, CRight),
   // the additional columns for (CRight, RExpand), and the diagonal (RExpand, RExpand).
   // We can add the rows first, and then add the columns (including the diagonal)
   // Alternatively, we can simply reconsruct the new F matrices from RNew
   StateComponent FRows = contract_from_right(herm(*HRight), RExpand, HamMatrices.right(), herm(R));
   F = tensor_col_sum(F, FRows);
   StateComponent FCols = contract_from_right(herm(*HRight), RNew, HamMatrices.right(), herm(RExpand));
   F = tensor_row_sum(F, FCols);
   HamMatrices.push_right(F);

   #else
   StateComponent RNew = tensor_col_sum(R, RExpand);
   OrthogonalizeBasis1_LQ(RNew);

   // reconstruct the F matrix with the expanded environment
   // if we reconstructed the complete F matrix, it would look like:
   HamMatrices.push_right(contract_from_right(herm(*HRight), RNew, HamMatrices.right(), herm(RNew)));

   // Update C in the new basis
   MatrixOperator U = scalar_prod(R, herm(RNew));
   *C = CC*U;
   *CRight = RNew;

   #endif

   return C->Basis2().total_dimension();
}

TruncationInfo DMRG::TruncateAndShiftLeft(StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector)
{
   TruncationInfo Info;
   MatrixOperator X = TruncateExpandBasis1(*C, HamMatrices.left(), *H, HamMatrices.right(), PostExpansionAlgo, States, ExtraStates, ExtraStatesPerSector, Info, Oversampling);
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
   MatrixOperator X = TruncateExpandBasis2(CC, HamMatrices.left(), *H, HamMatrices.right(), PostExpansionAlgo, States, ExtraStates, ExtraStatesPerSector, Info, Oversampling);
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

std::pair<MatrixOperator, RealDiagonalOperator>
SubspaceExpandBasis1(StateComponent& C, OperatorComponent const& H, StateComponent const& RightHam,
                     double MixFactor, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& LeftHam)
{
#if defined(SSC)
   MatrixOperator Lambda;
   SimpleStateComponent CX;
   std::tie(Lambda, CX) = ExpandBasis1_(C);
#else
   MatrixOperator Lambda = ExpandBasis1(C);
#endif

   MatrixOperator Rho = scalar_prod(herm(Lambda), Lambda);
   if (MixFactor > 0)
   {
#if defined(SSC)
      StateComponent RH = contract_from_right(herm(H), CX, RightHam, herm(CX));
#else
      StateComponent RH = contract_from_right(herm(H), C, RightHam, herm(C));
#endif
      MatrixOperator RhoMix;
      MatrixOperator RhoL = scalar_prod(Lambda, herm(Lambda));

      // Skip the identity and the Hamiltonian
      for (unsigned i = 1; i < RH.size()-1; ++i)
      {
         double Prefactor = norm_frob_sq(herm(Lambda) * LeftHam[i]);
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(herm(RH[i]), Rho, RH[i]);
	 //TRACE(i)(Prefactor);
      }
      // check for a zero mixing term - can happen if there are no interactions that span
      // the current bond
      double RhoTrace = trace(RhoMix).real();
      if (RhoTrace != 0)
         Rho += (MixFactor / RhoTrace) * RhoMix;
   }

   DensityMatrix<MatrixOperator> DM(Rho);
   //DM.DensityMatrixReport(std::cout);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);

   std::list<EigenInfo> KeptStates(DM.begin(), DMPivot);
   std::list<EigenInfo> DiscardStates(DMPivot, DM.end());

   MatrixOperator UKeep = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());
   Lambda = Lambda * herm(UKeep);

   //TRACE(Lambda);

#if defined(SSC)
   C = UKeep*CX; //prod(U, CX);
#else
   C = prod(UKeep, C);
#endif

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SingularValueDecompositionKeepBasis2(Lambda, U, D, Vh);

   //TRACE(U)(D)(Vh);

   C = prod(Vh, C);
   return std::make_pair(U, D);
}

// Apply subspace expansion / truncation on the right (C.Basis2()).
// Returns Lambda matrix (diagonal) and a unitary matrix
// Postcondition: C' Lambda' U' = C (up to truncation!)
std::pair<RealDiagonalOperator, MatrixOperator>
SubspaceExpandBasis2(StateComponent& C, OperatorComponent const& H, StateComponent const& LeftHam,
                     double MixFactor, StatesInfo const& States, TruncationInfo& Info,
                     StateComponent const& RightHam)
{
   MatrixOperator Lambda = ExpandBasis2(C);

   MatrixOperator Rho = scalar_prod(Lambda, herm(Lambda));
   if (MixFactor > 0)
   {
      StateComponent LH = contract_from_left(H, herm(C), LeftHam, C);
      MatrixOperator RhoMix;

      MatrixOperator RhoR = scalar_prod(herm(Lambda), Lambda);

      for (unsigned i = 1; i < LH.size()-1; ++i)
      {
         double Prefactor = norm_frob_sq(Lambda * RightHam[i]);
         if (Prefactor == 0)
            Prefactor = 1;
         RhoMix += Prefactor * triple_prod(LH[i], Rho, herm(LH[i]));
	 //	 TRACE(i)(Prefactor);
      }
      double RhoTrace = trace(RhoMix).real();
      if (RhoTrace != 0)
         Rho += (MixFactor / RhoTrace) * RhoMix;
   }
   DensityMatrix<MatrixOperator> DM(Rho);
   DensityMatrix<MatrixOperator>::const_iterator DMPivot =
      TruncateFixTruncationErrorRelative(DM.begin(), DM.end(),
                                         States,
                                         Info);
   std::list<EigenInfo> KeptStates(DM.begin(), DMPivot);
   std::list<EigenInfo> DiscardStates(DMPivot, DM.end());

   MatrixOperator UKeep = DM.ConstructTruncator(KeptStates.begin(), KeptStates.end());

   Lambda = UKeep * Lambda;
   C = prod(C, herm(UKeep));

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SingularValueDecompositionKeepBasis1(Lambda, U, D, Vh);

   C = prod(C, U);

   return std::make_pair(D, Vh);
}

TruncationInfo DMRG::TruncateAndShiftLeft3S(StatesInfo const& States, double MixFactor)
{
   MatrixOperator U;
   RealDiagonalOperator Lambda;
   TruncationInfo Info;
   LinearWavefunction::const_iterator CNext = C;
   --CNext;
   std::tie(U, Lambda) = SubspaceExpandBasis1(*C, *H, HamMatrices.right(), MixFactor, States, Info, HamMatrices.left());

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

   *C = prod(*C, U*Lambda);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

TruncationInfo DMRG::TruncateAndShiftRight3S(StatesInfo const& States, double MixFactor)
{
   // Truncate right
   RealDiagonalOperator Lambda;
   MatrixOperator U;
   TruncationInfo Info;
   LinearWavefunction::const_iterator CNext = C;
   ++CNext;
   std::tie(Lambda, U) = SubspaceExpandBasis2(*C, *H, HamMatrices.left(), MixFactor, States, Info, HamMatrices.right());
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

   *C = prod(Lambda*U, *C);

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
