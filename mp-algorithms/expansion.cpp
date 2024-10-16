// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/expansion.cpp
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
#include <fstream>

double const PrefactorEpsilon = 1e-16;

using MessageLogger::Logger;
using MessageLogger::msg_log;

// These are not needed in C++17 (?), but C++14 requires the initializer
constexpr std::array<char const*,5> PreExpansionTraits::Names;
constexpr std::array<char const*, 6> PostExpansionTraits::Names;

std::ofstream ExpandFile(getenv_or_default("MP_EXPANDFILE", ""), std::ios_base::out | std::ios_base::trunc);

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
   double Max = 0.0;
   for (auto m : w)
   {
      if (m > 0.0 && (m < Min || Min == 0.0))
         Min = m;
      if (m > Max)
         Max = m;
   }
   // If all of the weights are identically zero, then set them to be the same
   if (Min == 0.0)
   {
      Min = 1.0;
      Max = 1.0;
   }
   for (auto& m : w)
   {
      if (m == 0.0)
         m = Min;
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

// Get a basis for use in the randomized SVD, selecting ExtraStates + Oversampling states out of the available states.
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
      CMatSVD DM(CMat, States.KeepZeroEigenvalues() ? CMatSVD::Right : CMatSVD::Both, CMatSVD::Right);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator VKeep = DM.ConstructRightVectors(DM.begin(), DMPivot);
      R = ReshapeFromBasis2(VKeep, C.LocalBasis(), C.Basis2());
      Info.ExtraStates_ = 0;
   }
   else if (Algo == PostExpansionAlgorithm::RSVD || Algo == PostExpansionAlgorithm::RangeFinding || Algo == PostExpansionAlgorithm::Mixing)
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

      if (Algo == PostExpansionAlgorithm::RangeFinding)
      {
         // For rangefinding, we already have the basis, just add it to R and orthogonalize it
         Info.ExtraStates_ = RNull.Basis1().total_dimension();
         R = RegularizeBasis1(tensor_col_sum(R, RNull));
         OrthogonalizeBasis1_LQ(R);  // additional orthogonalization step, ensure vectors are all orthogonal
      }
      else if (Algo == PostExpansionAlgorithm::Mixing)
      {
         // The 'mixing' strategy:
         // First round SVD to m+k states.
         // Construct another k states in the tangent space of the first m states (or tangent space of m+k?)
         // Mix these into the first round matrix with a weight equal to g * weight of final k states, for some constant g
         double Trunc = Info.TruncationError();
         if (Trunc == 0)
            Trunc = Info.SmallestKeptEigenvalue();
         //TRACE(Trunc);
         // In the mixing strategy, we construct the SVD of the combined CMat with the 3S term
         StateComponent D = contract_from_right(herm(W*P*H), C, RightHam, herm(RNull));
         MatrixOperator X = ReshapeBasis1(D);
         // X is wdm × k
         // Need to do a factorization to make it k × k
         X = QR_Factorize(X).second;
         StateComponent XNull = X * RNull;
         X = ReshapeBasis2(XNull);
         double g = 1.0;
         X *= g*std::sqrt(Trunc)/norm_frob(X);
         MatrixOperator CX = tensor_col_sum(CMat, X);
         CMatSVD ExpandDM(CX, CMatSVD::Right);
         auto DMPivot = TruncateFixTruncationErrorRelative(ExpandDM.begin(), ExpandDM.end(), States, Info);
         auto AdditionalStates = TruncateExtraStates(DMPivot, ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
         Info.ExtraStates_ = AdditionalStates.size();
         // make the list of m+k states
         std::vector<EigenInfo> KeptStates(ExpandDM.begin(), DMPivot);
         KeptStates.insert(KeptStates.end(), AdditionalStates.begin(), AdditionalStates.end());
         // Construct the new kept states
         MatrixOperator QExpand = ExpandDM.ConstructRightVectors(KeptStates.begin(), KeptStates.end());
         R = ReshapeFromBasis2(QExpand, C.LocalBasis(), C.Basis2());
      }
      else if (Algo == PostExpansionAlgorithm::RSVD)
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

         Info.ExtraStates_ = RNull.Basis1().total_dimension();
         R = RegularizeBasis1(tensor_col_sum(R, RNull));
         OrthogonalizeBasis1_LQ(R);  // additional orthogonalization step, ensure vectors are all orthogonal
      }
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
      CMatSVD DM(CMat, States.KeepZeroEigenvalues() ? CMatSVD::Left : CMatSVD::Both, CMatSVD::Left);
      auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
      MatrixOperator UKeep = DM.ConstructLeftVectors(DM.begin(), DMPivot);
      L = ReshapeFromBasis1(UKeep, C.LocalBasis(), C.Basis1());
      Info.ExtraStates_ = 0;
   }
   else if (Algo == PostExpansionAlgorithm::RSVD || Algo == PostExpansionAlgorithm::RangeFinding || Algo == PostExpansionAlgorithm::Mixing)
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

      if (Algo == PostExpansionAlgorithm::RangeFinding)
      {
         Info.ExtraStates_ = LNull.Basis2().total_dimension();
         // LNull might not be exactly orthogonal
         L = RegularizeBasis2(tensor_row_sum(L, LNull));
         OrthogonalizeBasis2_QR(L);
      }
      else if (Algo == PostExpansionAlgorithm::Mixing)
      {
         // The 'mixing' strategy:
         // First round SVD to m+k states.
         // Construct another k states in the tangent space of the first m states (or tangent space of m+k?)
         // Mix these into the first round matrix with a weight equal to g * weight of final k states, for some constant g
         double Trunc = Info.TruncationError();
         if (Trunc == 0)
            Trunc = Info.SmallestKeptEigenvalue();
         // In the mixing strategy, we construct the SVD of the combined CMat with the 3S term
         StateComponent D = contract_from_left(H*herm(P)*W, herm(LNull), LeftHam, C);
         MatrixOperator X = ReshapeBasis2(D);
         X = LQ_Factorize(X).first;
         StateComponent XNull = LNull * X;
         X = ReshapeBasis1(XNull);
         double g = 1.0;
         X *= g*std::sqrt(Trunc)/norm_frob(X);
         MatrixOperator CX = tensor_row_sum(CMat, X);
         CMatSVD ExpandDM(CX, CMatSVD::Left);
         auto DMPivot = TruncateFixTruncationErrorRelative(ExpandDM.begin(), ExpandDM.end(), States, Info);
         auto AdditionalStates = TruncateExtraStates(DMPivot, ExpandDM.end(), ExtraStates, ExtraStatesPerSector, false);
         Info.ExtraStates_ = AdditionalStates.size();
         // make the list of m+k states
         std::vector<EigenInfo> KeptStates(ExpandDM.begin(), DMPivot);
         KeptStates.insert(KeptStates.end(), AdditionalStates.begin(), AdditionalStates.end());
         // Construct the new kept states
         MatrixOperator QExpand = ExpandDM.ConstructLeftVectors(KeptStates.begin(), KeptStates.end());
         L = ReshapeFromBasis1(QExpand, C.LocalBasis(), C.Basis1());
      }
      else if (Algo == PostExpansionAlgorithm::RSVD)
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
         Info.ExtraStates_ = LNull.Basis2().total_dimension();
         // LNull might not be exactly orthogonal
         L = RegularizeBasis2(tensor_row_sum(L, LNull));
         OrthogonalizeBasis2_QR(L);
      }
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

      if (ExpandFile)
      {
         ExpandFile << '\n';
         LinearAlgebra::Vector<double> v = SingularValues(X);
         for (int i = 0; i < v.size(); ++i)
         {
            ExpandFile << v[i] << '\n';
         }
      }

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
