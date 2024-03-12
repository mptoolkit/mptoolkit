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

DMRG::DMRG(int Verbose_)
   : ProjectTwoSiteTangent(false),
   TotalNumSweeps(0),
   TotalNumIterations(0),
   TotalNumMultiplies(0),
   LastSweepTime(0),
   LastSweepFidelity(0),
   Verbose(Verbose_)
{
}

void
DMRG::InitializeLeftOrtho(LinearWavefunction Psi_, BasicTriangularMPO const& Ham_, StateComponent const& E, StateComponent const& F)
{
   Hamiltonian = Ham_;
   Psi = std::move(Psi_);

   // construct the HamMatrix elements
   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   // Initialize the E, F matrices

   HamMatrices.push_left(E);

   H = Hamiltonian.begin();
   C = Psi.begin();
   Site_ = 0;
   for (int i = 0; i < Psi.size()-1; ++i)
   {
      if (Verbose > 1)
         std::cout << "site " << (HamMatrices.size_left()) << std::endl;
      HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));
      ++H;
      ++C;
      ++Site_;
   }

   SweepC = *C;

   HamMatrices.push_right(F);

   // clear the global statistics
   TotalNumSweeps = 0;
   TotalNumIterations = 0;
   TotalNumMultiplies = 0;

  this->debug_check_structure();
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
   ++SweepNumIterations;
   SweepNumMultiplies += IterationNumMultiplies;
   SweepSumStates += IterationNumStates;
   SweepMaxStates = std::max(SweepMaxStates, IterationNumStates);
   SweepMinEnergy = std::min(SweepMinEnergy, IterationEnergy.real());
   SweepTotalTruncation += IterationTruncation;
   SweepMaxEntropy = std::max(SweepMaxEntropy, IterationEntropy);
   SweepEnergyEachIteration.push_back(IterationEnergy);
   SweepEnergyVariance = statistics::variance(SweepEnergyEachIteration.begin(), SweepEnergyEachIteration.end());
}

void DMRG::StartSweep(bool IncrementSweepNumber, double /* Broad_ */)
{
   ++TotalNumSweeps;

   SweepNumIterations = 0;
   SweepSumStates = 0;
   SweepMaxStates = 0;
   SweepNumMultiplies = 0;
   SweepMinEnergy = 1e100;
   SweepTotalTruncation = 0;
   SweepMaxEntropy = 0;
   SweepStartTime = ProcControl::GetCumulativeElapsedTime();
   SweepEnergyVariance = std::numeric_limits<double>::infinity();
   SweepEnergyEachIteration.clear();

   SweepC = *C;
}

void DMRG::EndSweep()
{
   TotalNumIterations += SweepNumIterations;
   TotalNumMultiplies += SweepNumMultiplies;

   LastSweepTime = ProcControl::GetCumulativeElapsedTime() - SweepStartTime;
   LastSweepFidelity = norm_frob(inner_prod(*C, SweepC));
}

std::complex<double>
DMRG::Energy() const
{
   return inner_prod(HamMatrices.right(), contract_from_left(*H, herm(HamMatrices.left()), *C, HamMatrices.right()));
}

std::complex<double>
DMRG::Solve()
{
   Solver_.Solve(*C, HamMatrices.left(), *H, HamMatrices.right());

   // NOTE: *C is typically nearly, but not exactly normalized at this point.
   // There is no point normalizing here, since we need to do it anyway after the truncation.

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   return IterationEnergy;
}

std::pair<std::complex<double>, TruncationInfo>
DMRG::SolveCoarseGrainRight(StatesInfo const& States)
{
   auto R = C; ++R;
   StateComponent Rold = *R;
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

   AMatSVD DM(X, CoarseGrainBasis);

   TruncationInfo Info;
   auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
   RealDiagonalOperator Lambda;
   DM.ConstructMatrices(DM.begin(), DMPivot, *C, Lambda, *R);

   Lambda *= 1.0 / norm_frob(Lambda); // normalize

   // update block Hamiltonian
   HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));

   MatrixOperator SweepCProj = scalar_prod(herm(*C), SweepC);

   // next site
   ++Site_;      // note: this is now one site ahead of the site that was optimized.
   ++C;
   SweepC = prod(SweepCProj, Rold);
   *C = Lambda * (*C);
   ++H;

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = Info.KeptEntropy();

   return std::make_pair(IterationEnergy, Info);
}

std::pair<std::complex<double>, TruncationInfo>
DMRG::SolveCoarseGrainLeft(StatesInfo const& States)
{
   auto L = C; --L;
   StateComponent Lold = *L;
   auto CoarseGrainBasis = make_product_basis(L->LocalBasis(), C->LocalBasis());
   StateComponent X = local_tensor_prod(*L, *C);

   auto HL = H; --HL;
   OperatorComponent HX = local_tensor_prod(*HL, *H);

   HamMatrices.pop_left();

   Solver_.Solve(X, HamMatrices.left(), HX, HamMatrices.right());

   IterationEnergy = Solver_.LastEnergy();
   IterationNumMultiplies = Solver_.LastIter();

   AMatSVD DM(X, CoarseGrainBasis);

   TruncationInfo Info;
   auto DMPivot = TruncateFixTruncationErrorRelative(DM.begin(), DM.end(), States, Info);
   RealDiagonalOperator Lambda;
   DM.ConstructMatrices(DM.begin(), DMPivot, *L, Lambda, *C);

   Lambda *= 1.0 / norm_frob(Lambda); // normalize

   // update block Hamiltonian
   HamMatrices.push_right(contract_from_right(herm(*H), *C, HamMatrices.right(), herm(*C)));

   MatrixOperator SweepCProj = scalar_prod(SweepC, herm(*C));

   // next site
   --Site_;         // note: this is now one site ahead of the site that was optimized.
   --C;             // Now C == L
   --H;
   SweepC = prod(Lold, SweepCProj);
   *C = (*C) * Lambda;

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = Info.KeptEntropy();

   return std::make_pair(IterationEnergy, Info);
}

int DMRG::ExpandLeftEnvironment(int StatesWanted, int ExtraStatesPerSector)
{
   // Quick return if there is no need to expand the basis
   if (PreExpansionAlgo == PreExpansionAlgorithm::NoExpansion || (StatesWanted <= C->Basis1().total_dimension() && ExtraStatesPerSector <= 0))
      return C->Basis1().total_dimension();

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
   SweepC = tensor_col_sum(SweepC, CExpand);

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
   SweepC = U * SweepC;

   #endif

   return C->Basis1().total_dimension();
}

int DMRG::ExpandRightEnvironment(int StatesWanted, int ExtraStatesPerSector)
{
   // Quick return if there is no need to expand the basis
   if (PreExpansionAlgo == PreExpansionAlgorithm::NoExpansion || (StatesWanted <= C->Basis2().total_dimension() && ExtraStatesPerSector <= 0))
      return C->Basis2().total_dimension();

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
   SweepC = tensor_row_sum(SweepC, CExpand);

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
   SweepC = SweepC*U;

   #endif

   return C->Basis2().total_dimension();
}

void
DMRG::ShiftRight(MatrixOperator const& Lambda)
{
   // update blocks
   HamMatrices.push_left(contract_from_left(*H, herm(*C), HamMatrices.left(), *C));
   HamMatrices.pop_right();

   // The projection operator that maps from the SweepC basis to the new basis
   MatrixOperator SweepCProj = scalar_prod(herm(*C), SweepC);

   // next site
   ++Site_;
   ++H;
   ++C;

   SweepC = prod(SweepCProj, *C);

   *C = prod(Lambda, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

}

void
DMRG::ShiftLeft(MatrixOperator const& Lambda)
{
   // update blocks
   HamMatrices.push_right(contract_from_right(herm(*H), *C, HamMatrices.right(), herm(*C)));
   HamMatrices.pop_left();

   // The projection operator that maps from the SweepC basis to the new basis
   MatrixOperator SweepCProj = scalar_prod(SweepC, herm(*C));

   // Move to the next site
   --Site_;
   --H;
   --C;

   SweepC = prod(*C, SweepCProj);

   *C = prod(*C, Lambda);

   // normalize
   *C *= 1.0 / norm_frob(*C);
}

TruncationInfo DMRG::TruncateAndShiftLeft(StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector)
{
   TruncationInfo Info;
   MatrixOperator X = TruncateExpandBasis1(*C, HamMatrices.left(), *H, HamMatrices.right(), PostExpansionAlgo, States, ExtraStates, ExtraStatesPerSector, Info, Oversampling);
   if (Verbose > 1)
   {
      std::cerr << "Truncating left basis, states=" << Info.KeptStates() << '\n';
   }

   this->ShiftLeft(X);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

TruncationInfo DMRG::TruncateAndShiftRight(StatesInfo const& States, int ExtraStates, int ExtraStatesPerSector)
{
   TruncationInfo Info;
   MatrixOperator X = TruncateExpandBasis2(*C, HamMatrices.left(), *H, HamMatrices.right(), PostExpansionAlgo, States, ExtraStates, ExtraStatesPerSector, Info, Oversampling);
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }

   this->ShiftRight(X);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

TruncationInfo DMRG::TruncateAndShiftLeft3S(StatesInfo const& States, double MixFactor)
{
   MatrixOperator U;
   RealDiagonalOperator Lambda;
   TruncationInfo Info;
   std::tie(U, Lambda) = SubspaceExpandBasis1(*C, *H, HamMatrices.right(), MixFactor, States, Info, HamMatrices.left());

   if (Verbose > 1)
   {
      std::cerr << "Truncating left basis, states=" << Info.KeptStates() << '\n';
   }

   this->ShiftLeft(U*Lambda);

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
   std::tie(Lambda, U) = SubspaceExpandBasis2(*C, *H, HamMatrices.left(), MixFactor, States, Info, HamMatrices.right());
   if (Verbose > 1)
   {
      std::cerr << "Truncating right basis, states=" << Info.KeptStates() << '\n';
   }

   this->ShiftRight(Lambda*U);

   IterationNumStates = Info.KeptStates();
   IterationTruncation += Info.TruncationError();
   IterationEntropy = std::max(IterationEntropy, Info.KeptEntropy());

   return Info;
}

void DMRG::check_structure() const
{
   CHECK_EQUAL(HamMatrices.left().Basis1(), C->Basis1());
   CHECK_EQUAL(HamMatrices.right().Basis1(), C->Basis2());

   CHECK_EQUAL(HamMatrices.left().Basis2(), C->Basis1());
   CHECK_EQUAL(HamMatrices.right().Basis2(), C->Basis2());

   CHECK_EQUAL(C->LocalBasis(), H->LocalBasis1());
   CHECK_EQUAL(H->LocalBasis1(), H->LocalBasis1());

   CHECK_EQUAL(C->LocalBasis(), SweepC.LocalBasis());
   CHECK_EQUAL(C->Basis1(), SweepC.Basis1());
   CHECK_EQUAL(C->Basis2(), SweepC.Basis2());

}
void DMRG::debug_check_structure() const
{
   this->check_structure();
   #if !defined(NDEBUG)
   #endif
}
