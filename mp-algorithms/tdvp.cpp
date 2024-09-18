// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/tdvp.cpp
//
// Copyright (C) 2021-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
// Copyright (C) 2023 Ian McCulloch <ian@qusim.net>
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

#include "tdvp.h"
#include "lanczos-exponential-new.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"
#include "lattice/infinite-parser.h"
#include "parser/parser.h"

double const PrefactorEpsilon = 1e-16;

std::complex<double> const I(0.0, 1.0);

struct HEff1
{
   HEff1(StateComponent const& E_, OperatorComponent const& H_, StateComponent const& F_)
      : E(E_), H(H_), F(F_)
   {
   }

   StateComponent operator()(StateComponent const& x) const
   {
      StateComponent R = operator_prod_inner(H, E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   OperatorComponent const& H;
   StateComponent const& F;
};

struct HEff2
{
   HEff2(StateComponent const& E_, StateComponent const& F_)
      : E(E_), F(F_)
   {
   }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator R = operator_prod(E, x, herm(F));
      return R;
   }

   StateComponent const& E;
   StateComponent const& F;
};

Hamiltonian::Hamiltonian(std::string HamStr, int Size_,
                         std::string Magnus_, std::string TimeVar_, int Verbose_)
   : Size(Size_), Magnus(Magnus_), TimeVar(TimeVar_), Verbose(Verbose_)
{
   std::tie(HamOperator, Lattice) = ParseOperatorStringAndLattice(HamStr);

   // Attempt to convert the operator into an MPO. If it works, then the
   // Hamiltonian is time-independent. If it fails, assume it failed because
   // the time variable isn't defined yet.  If there is some other error in the
   // operator, we'll catch it later.
   try
   {
      HamMPO = ParseTriangularOperator(Lattice, HamOperator);
      if (Size != 0)
         HamMPO = repeat(HamMPO, Size / HamMPO.size());
      TimeDependent = false;
   }
   catch (Parser::ParserError& e)
   {
      if (Verbose > 1)
         std::cerr << "Parser error converting the Hamiltonian to an MPO - assuming the Hamiltonian is time-dependent." << std::endl;
      TimeDependent = true;
   }
   catch (...)
   {
      throw;
   }
}

BasicTriangularMPO
Hamiltonian::operator()(std::complex<double> t, std::complex<double> dt) const
{
   if (TimeDependent == false)
      return HamMPO;
   else
   {
      BasicTriangularMPO HamMPO_;

      if (dt == 0.0)
         HamMPO_ = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, t}});
      else if (Magnus == "2")
      {
         std::complex<double> Time = t + 0.5*dt;
         HamMPO_ = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, Time}});
      }
      else if (Magnus == "4")
      {
         std::complex<double> t1 = t + (0.5 - std::sqrt(3.0)/6.0) * dt;
         std::complex<double> t2 = t + (0.5 + std::sqrt(3.0)/6.0) * dt;
         BasicTriangularMPO A1 = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, t1}});
         BasicTriangularMPO A2 = ParseTriangularOperator(Lattice, HamOperator, {{TimeVar, t2}});
         HamMPO_ = 0.5 * (A1 + A2) - I * dt * (std::sqrt(3.0) / 12.0) * commutator(A1, A2);
      }

      if (Size != 0)
         HamMPO_ = repeat(HamMPO_, Size / HamMPO_.size());

      return HamMPO_;
   }
}

void
Hamiltonian::set_size(int Size_)
{
   Size = Size_;

   if (Size != 0 && TimeDependent == false)
      HamMPO = repeat(HamMPO, Size / HamMPO.size());
}

TDVP::TDVP(Hamiltonian const& Ham_, TDVPSettings const& Settings_)
   : Ham(Ham_), InitialTime(Settings_.InitialTime), Timestep(Settings_.Timestep),
     Comp(Settings_.Comp), MaxIter(Settings_.MaxIter), ErrTol(Settings_.ErrTol),
     SInfo(Settings_.SInfo), ExpandFactor(Settings_.ExpandFactor),
     ExpandMinStates(Settings_.ExpandMinStates), ExpandMinPerSector(Settings_.ExpandMinPerSector),
     Epsilon(Settings_.Epsilon), Normalize(Settings_.Normalize), Verbose(Settings_.Verbose)
{
}

TDVP::TDVP(FiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_, TDVPSettings const& Settings_)
   : TDVP(Ham_, Settings_)
{
   // Initialize Psi and Ham.
   Time = InitialTime;
   std::complex<double> dt = Comp.Beta.back()*Timestep;
   HamMPO = Ham(Time-dt, dt);

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;
   H = HamMPO.begin();
   HamL.push_back(Initial_E(HamMPO, Psi_.Basis1()));
   for (FiniteWavefunctionLeft::const_mps_iterator I = Psi_.begin(); I != Psi_.end(); ++I)
   {
      if (Verbose > 1)
         std::cout << "Site " << (HamL.size()) << std::endl;
      HamL.push_back(contract_from_left(*H, herm(*I), HamL.back(), *I));
      Psi.push_back(*I);
      MaxStates = std::max(MaxStates, I->Basis2().total_dimension());
      ++H;
   }
   HamR.push_front(Initial_F(HamMPO, Psi_.Basis2()));
   Psi.set_back(prod(Psi.get_back(), Psi_.lambda_r()));

   // Initialize to the right-most site.
   HamL.pop_back();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size() - 1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;
}

FiniteWavefunctionLeft
TDVP::Wavefunction() const
{
   return FiniteWavefunctionLeft::Construct(Psi);
}

std::complex<double>
TDVP::Energy() const
{
   return inner_prod(contract_from_left(*H, herm(*C), HamL.back(), *C), HamR.front());
}

void
TDVP::EvolveCurrentSite(std::complex<double> Tau)
{
   // Evolve current site.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, -I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void
TDVP::IterateLeft(std::complex<double> Tau)
{
   // Perform SVD to right-orthogonalize current site.
   MatrixOperator U;
   RealDiagonalOperator D;

   std::tie(U, D) = OrthogonalizeBasis1(*C);

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

   // Evolve the UD term backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;
   MatrixOperator UD = U*D;

   UD = LanczosExponential(UD, HEff2(HamL.back(), HamR.front()), Iter, I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " UD Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Move to the next site.
   --Site;
   --H;
   --C;

   *C = prod(*C, UD);

   HamL.pop_back();
}

void
TDVP::IterateRight(std::complex<double> Tau)
{
   // Perform SVD to left-orthogonalize current site.
   MatrixOperator Vh;
   RealDiagonalOperator D;

   std::tie(D, Vh) = OrthogonalizeBasis2(*C);

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));

   // Evolve the DVh term backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;
   MatrixOperator DVh = D*Vh;

   DVh = LanczosExponential(DVh, HEff2(HamL.back(), HamR.front()), Iter, I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " DVh Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Move to the next site.
   ++Site;
   ++H;
   ++C;

   *C = prod(DVh, *C);

   HamR.pop_front();
}

// ***
// Borrowed from mp-algorithms.dmrg
// ***

// Calculate how to determine ExtraStates among the available states, with relative weights given by the Weights array.
// To do this, we start from the total weight of the available sectors, and stretch that interval w[j] to be length ExtraStates.
// (we can cap ExtraStates at the total number of available states in sectors with non-zero weight, so ExtraStates is actually
// achievable. We can also remove sectors that have zero weight from the Avail array.)
// Then the number of states we keep of the j'th sector, k[j] is k[j] = round(sum(k<=j) w[k] - k[j-1]) (with k[-1] = 0).
// If any of the k[j] > Avail[j], then pin k[j] = Avail[j], remove sector j from the Avail array, and repeat the distribution.
// We avoid creating subspaces that have zero dimension. In principle this is harmless, although inefficient, but
// MKL doesn't behave properly when multiplying 0-size matrices.
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
         // randomize the order of the sectors
         randutil::random_stream stream;
         stream.seed();
         std::shuffle(ToAdd.begin(), ToAdd.end(), stream.u_rand);
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

void
ExpandLeftEnvironment(StateComponent& CLeft, StateComponent& CRight,
                      StateComponent const& E, StateComponent const& F,
                      OperatorComponent const& HLeft, OperatorComponent const& HRight,
                      StatesInfo SInfo, double ExpandFactor, int ExpandMinStates,
                      int ExpandMinPerSector, int Verbose)
{
   // If the MPO dimension across the bond is two, we cannot do any expansion anyway.
   if (HLeft.Basis2().size() == 2)
      return;

   // Perform truncation before expansion.
   CMatSVD SVD(ExpandBasis1(CRight));
   TruncationInfo Info;
   auto Cutoff = TruncateFixTruncationErrorRelative(SVD.begin(), SVD.end(), SInfo, Info);

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SVD.ConstructMatrices(SVD.begin(), Cutoff, U, D, Vh);

   if (Verbose > 0)
      std::cout << "StatesOld=" << CLeft.Basis2().total_dimension()
                << " StatesTrunc=" << Info.KeptStates();

   CRight = prod(D*Vh, CRight);
   CLeft = prod(CLeft, U);

   // Calculate left null space of left site.
   StateComponent NLeft = NullSpace2(CLeft);

   // Perform SVD to right-orthogonalize the right site and extract singular value matrix.
   StateComponent CRightOrtho = CRight;
   MatrixOperator URight;
   RealDiagonalOperator DRight;
   std::tie(URight, DRight) = OrthogonalizeBasis1(CRightOrtho);

   // Project out the first and last columns of HLeft.
   SimpleOperator Projector(HLeft.Basis2(), BasisList(HLeft.Basis2().begin()+1, HLeft.Basis2().end()-1));
   for (int i = 0; i < Projector.Basis2().size(); ++i)
      Projector(i+1, i) = 1.0;

   StateComponent X = contract_from_left(HLeft*Projector, herm(NLeft), E, CLeft*URight*DRight);
   StateComponent FRight = contract_from_right(herm(herm(Projector)*HRight), CRightOrtho, F, herm(CRight));

   // Multiply each element of X by a prefactor depending on the corresponding element of F.
   for (int i = 0; i < X.size(); ++i)
   {
      // We add an epsilon term to the prefactor to avoid the corner case where
      // the prefactor is zero for all elements, and any possible new states
      // are ignored.
      double Prefactor = norm_frob(FRight[i]) + PrefactorEpsilon;
      X[i] *= Prefactor;
   }

   MatrixOperator XExpand = ExpandBasis1(X);

   // Randomized SVD.
   std::map<QuantumNumbers::QuantumNumber, int> NumAvailablePerSector = RankPerSector(XExpand);
   std::map<QuantumNumbers::QuantumNumber, int> KeptStateDimension = DimensionPerSector(CLeft.Basis2());
   std::map<QuantumNumbers::QuantumNumber, int> Weights;
   for (auto n : NumAvailablePerSector)
   {
      // Set the relative weight to the number of kept states in this sector.
      // Clamp this at a minimum of 1 state, since a state only appears in X if it has non-zero contribution
      // with respect to the density matrix mixing.
      if (n.second > 0)
         Weights[n.first] = std::max(KeptStateDimension[n.first], 1);
   }

   // Target number of extra states.
   int ExtraStates = std::max((int) std::ceil(ExpandFactor * CLeft.Basis2().total_dimension()), ExpandMinStates);

   auto NumExtraPerSector = DistributeStates(Weights, NumAvailablePerSector, ExtraStates*2, ExpandMinPerSector*2);

   VectorBasis ExtraBasis(CLeft.GetSymmetryList(), NumExtraPerSector.begin(), NumExtraPerSector.end());

   MatrixOperator M = MakeRandomMatrixOperator(XExpand.Basis2(), ExtraBasis);
   MatrixOperator A = XExpand*M;
   auto QR = QR_Factorize(A);
   XExpand = herm(QR.first) * XExpand;

   // Take the truncated SVD of X.
   CMatSVD SVDEx(XExpand);
   auto ExpandedStates = TruncateExtraStates(SVDEx.begin(), SVDEx.end(), ExtraStates, ExpandMinPerSector, false);

   MatrixOperator UEx, VhEx;
   RealDiagonalOperator DEx;
   SVDEx.ConstructMatrices(ExpandedStates.begin(), ExpandedStates.end(), UEx, DEx, VhEx);

   // Construct new basis.
   SumBasis<VectorBasis> NewBasis(CLeft.Basis2(), UEx.Basis2());
   // Regularize the new basis.
   Regularizer R(NewBasis);

   // Add the new states to CLeft, and add zeros to CRight.
   CLeft = RegularizeBasis2(tensor_row_sum(CLeft, prod(NLeft, QR.first * UEx), NewBasis), R);

   StateComponent Z = StateComponent(CRight.LocalBasis(), VhEx.Basis1(), CRight.Basis2());
   CRight = RegularizeBasis1(R, tensor_col_sum(CRight, Z, NewBasis));

   if (Verbose > 0)
      std::cout << " StatesNew=" << CLeft.Basis2().total_dimension() << std::endl;
}

void
ExpandRightEnvironment(StateComponent& CLeft, StateComponent& CRight,
                       StateComponent const& E, StateComponent const& F,
                       OperatorComponent const& HLeft, OperatorComponent const& HRight,
                       StatesInfo SInfo, double ExpandFactor, int ExpandMinStates,
                       int ExpandMinPerSector, int Verbose)
{
   // If the MPO dimension across the bond is two, we cannot do any expansion anyway.
   if (HRight.Basis1().size() == 2)
      return;

   // Perform truncation before expansion.
   CMatSVD SVD(ExpandBasis2(CLeft));
   TruncationInfo Info;
   auto Cutoff = TruncateFixTruncationErrorRelative(SVD.begin(), SVD.end(), SInfo, Info);

   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SVD.ConstructMatrices(SVD.begin(), Cutoff, U, D, Vh);

   if (Verbose > 0)
      std::cout << "StatesOld=" << CRight.Basis1().total_dimension()
                << " StatesTrunc=" << Info.KeptStates();

   CRight = prod(Vh, CRight);
   CLeft = prod(CLeft, U*D);

   // Calculate right null space of right site.
   StateComponent NRight = NullSpace1(CRight);

   // Perform SVD to left-orthogonalize the left site and extract singular value matrix.
   StateComponent CLeftOrtho = CLeft;
   MatrixOperator VhLeft;
   RealDiagonalOperator DLeft;
   std::tie(DLeft, VhLeft) = OrthogonalizeBasis2(CLeftOrtho);

   // Project out the first and last columns of HRight.
   SimpleOperator Projector(BasisList(HRight.Basis1().begin()+1, HRight.Basis1().end()-1), HRight.Basis1());
   for (int i = 0; i < Projector.Basis1().size(); ++i)
      Projector(i, i+1) = 1.0;

   StateComponent X = contract_from_right(herm(Projector*HRight), NRight, F, herm(DLeft*VhLeft*CRight));
   StateComponent ELeft = contract_from_left(HLeft*herm(Projector), herm(CLeftOrtho), E, CLeft);

   // Multiply each element of X by a prefactor depending on the corresponding element of E.
   for (int i = 0; i < X.size(); ++i)
   {
      // We add an epsilon term to the prefactor to avoid the corner case where
      // the prefactor is zero for all elements, and any possible new states
      // are ignored.
      double Prefactor = norm_frob(ELeft[i]) + PrefactorEpsilon;
      X[i] *= Prefactor;
   }

   MatrixOperator XExpand = ExpandBasis1(X);

   // Randomized SVD.
   std::map<QuantumNumbers::QuantumNumber, int> NumAvailablePerSector = RankPerSector(XExpand);
   std::map<QuantumNumbers::QuantumNumber, int> KeptStateDimension = DimensionPerSector(CLeft.Basis2());
   std::map<QuantumNumbers::QuantumNumber, int> Weights;
   for (auto n : NumAvailablePerSector)
   {
      // Set the relative weight to the number of kept states in this sector.
      // Clamp this at a minimum of 1 state, since a state only appears in X if it has non-zero contribution
      // with respect to the density matrix mixing.
      if (n.second > 0)
         Weights[n.first] = std::max(KeptStateDimension[n.first], 1);
   }

   // Target number of extra states.
   int ExtraStates = std::max((int) std::ceil(ExpandFactor * CLeft.Basis2().total_dimension()), ExpandMinStates);

   auto NumExtraPerSector = DistributeStates(Weights, NumAvailablePerSector, ExtraStates*2, ExpandMinPerSector*2);

   VectorBasis ExtraBasis(CLeft.GetSymmetryList(), NumExtraPerSector.begin(), NumExtraPerSector.end());

   MatrixOperator M = MakeRandomMatrixOperator(XExpand.Basis2(), ExtraBasis);
   MatrixOperator A = XExpand*M;
   auto QR = QR_Factorize(A);
   XExpand = herm(QR.first) * XExpand;

   // Take the truncated SVD of X.
   CMatSVD SVDEx(XExpand);
   auto ExpandedStates = TruncateExtraStates(SVDEx.begin(), SVDEx.end(), ExtraStates, ExpandMinPerSector, false);

   MatrixOperator UEx, VhEx;
   RealDiagonalOperator DEx;
   SVDEx.ConstructMatrices(ExpandedStates.begin(), ExpandedStates.end(), UEx, DEx, VhEx);

   // Construct new basis.
   SumBasis<VectorBasis> NewBasis(CRight.Basis1(), UEx.Basis2());
   // Regularize the new basis.
   Regularizer R(NewBasis);

   // Add the new states to CRight, and add zeros to CLeft.
   CRight = RegularizeBasis1(R, tensor_col_sum(CRight, prod(herm(QR.first * UEx), NRight), NewBasis));

   StateComponent Z = StateComponent(CLeft.LocalBasis(), CLeft.Basis1(), VhEx.Basis1());
   CLeft = RegularizeBasis2(tensor_row_sum(CLeft, Z, NewBasis), R);

   if (Verbose > 0)
      std::cout << " StatesNew=" << CRight.Basis1().total_dimension() << std::endl;
}

void
TDVP::ExpandLeft()
{
   auto CNext = C;
   --CNext;
   auto HNext = H;
   --HNext;

   HamL.pop_back();

   if (Verbose > 1)
      std::cout << "Timestep=" << TStep
                << " Site=" << Site << " ";

   ExpandLeftEnvironment(*CNext, *C, HamL.back(), HamR.front(), *HNext, *H,
                         SInfo, ExpandFactor, ExpandMinStates, ExpandMinPerSector, Verbose-1);

   MaxStates = std::max(MaxStates, C->Basis1().total_dimension());

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*HNext, herm(*CNext), HamL.back(), *CNext));
}

void
TDVP::ExpandRight()
{
   auto CNext = C;
   ++CNext;
   auto HNext = H;
   ++HNext;

   HamR.pop_front();

   if (Verbose > 1)
      std::cout << "Timestep=" << TStep
                << " Site=" << Site << " ";

   ExpandRightEnvironment(*C, *CNext, HamL.back(), HamR.front(), *H, *HNext,
                          SInfo, ExpandFactor, ExpandMinStates, ExpandMinPerSector, Verbose-1);

   MaxStates = std::max(MaxStates, C->Basis2().total_dimension());

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*HNext), *CNext, HamR.front(), herm(*CNext)));
}

void
TDVP::SweepLeft(std::complex<double> Tau, bool Expand)
{
   while (Site > LeftStop)
   {
      if (Expand)
         this->ExpandLeft();
      this->EvolveCurrentSite(Tau);
      this->IterateLeft(Tau);
   }

   this->EvolveCurrentSite(Tau);
}

void
TDVP::SweepRight(std::complex<double> Tau, bool Expand)
{
   while (Site < RightStop)
   {
      if (Expand)
         this->ExpandRight();
      this->EvolveCurrentSite(Tau);
      this->IterateRight(Tau);
   }

   this->EvolveCurrentSite(Tau);
}

void
TDVP::CalculateEps1()
{
   // Perform SVD to right-orthogonalize current site for NullSpace1.
   StateComponent CRightOrtho = *C;
   OrthogonalizeBasis1(CRightOrtho);

   // Calculate error measure epsilon_1 and add to sum.
   StateComponent Y = contract_from_right(herm(*H), NullSpace1(CRightOrtho), HamR.front(), herm(*C));
   double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y)));
   Eps1SqSum += Eps1Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " Eps1Sq=" << Eps1Sq
                << std::endl;
   }
}

void
TDVP::CalculateEps12()
{
   // Calculate the left half of epsilon_2.
   StateComponent HamLBack = HamL.back();
   HamL.pop_back();

   LinearWavefunction::const_iterator CPrev = C;
   --CPrev;
   BasicTriangularMPO::const_iterator HPrev = H;
   --HPrev;

   StateComponent X = contract_from_left(*HPrev, herm(NullSpace2(*CPrev)), HamL.back(), *CPrev);

   HamL.push_back(HamLBack);

   // Perform SVD to right-orthogonalize current site for NullSpace1.
   StateComponent CRightOrtho = *C;
   OrthogonalizeBasis1(CRightOrtho);

   // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
   StateComponent Y = contract_from_right(herm(*H), NullSpace1(CRightOrtho), HamR.front(), herm(*C));
   double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y)));
   double Eps2Sq = norm_frob_sq(scalar_prod(X, herm(Y)));
   Eps1SqSum += Eps1Sq;
   Eps2SqSum += Eps2Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " Eps1Sq=" << Eps1Sq
                << " Eps2Sq=" << Eps2Sq
                << std::endl;
   }
}

void
TDVP::SweepRightFinal(std::complex<double> Tau, bool Expand)
{
   if (Expand && Site < RightStop)
      this->ExpandRight();
   this->EvolveCurrentSite(Tau);
   this->CalculateEps1();

   while (Site < RightStop)
   {
      this->IterateRight(Tau);
      if (Expand && Site < RightStop)
         this->ExpandRight();
      this->EvolveCurrentSite(Tau);
      this->CalculateEps12();
   }
}

void
TDVP::Evolve(bool Expand)
{
   // Reset MaxStates if we are expanding.
   if (Expand)
      MaxStates = 0;

   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Alpha = Comp.Alpha.cbegin();
   std::vector<double>::const_iterator Beta = Comp.Beta.cbegin();

   this->UpdateHamiltonianLeft(Time, (*Alpha)*Timestep);
   this->SweepLeft((*Alpha)*Timestep, Expand);
   Time += (*Alpha)*Timestep;
   ++Alpha;

   while (Alpha != Comp.Alpha.cend())
   {
      this->UpdateHamiltonianRight(Time, (*Beta)*Timestep);
      this->SweepRight((*Beta)*Timestep, Expand);
      Time += (*Beta)*Timestep;
      ++Beta;

      this->UpdateHamiltonianLeft(Time, (*Alpha)*Timestep);
      this->SweepLeft((*Alpha)*Timestep, Expand);
      Time += (*Alpha)*Timestep;
      ++Alpha;
   }

   this->UpdateHamiltonianRight(Time, (*Beta)*Timestep);

   if (Epsilon)
      this->SweepRightFinal((*Beta)*Timestep, Expand);
   else
      this->SweepRight((*Beta)*Timestep, Expand);

   Time += (*Beta)*Timestep;
}

void
TDVP::IterateLeft2(std::complex<double> Tau)
{
   // Form two-site centre block.
   --Site;

   LinearWavefunction::iterator CPrev = C;
   --C;
   StateComponent C2 = local_tensor_prod(*C, *CPrev);

   BasicTriangularMPO::const_iterator HPrev = H;
   --H;
   OperatorComponent H2 = local_tensor_prod(*H, *HPrev);

   HamL.pop_back();

   // Evolve two-site centre block.
   int Iter = MaxIter;
   double Err = ErrTol;

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, -I*Tau, Err, LogAmplitude);

   // Perform SVD on new C2.
   AMatSVD SL(C2, Tensor::ProductBasis<BasisList, BasisList>(C->LocalBasis(), CPrev->LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationErrorAbsolute(SL.begin(), SL.end(), SInfo, Info);
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, *C, D, *CPrev);
   *C = prod(*C, D);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Sites=" << Site << "," << Site+1
                << " C2 Iter=" << Iter
                << " Err=" << Err
                << " States=" << Info.KeptStates()
		<< " TruncErr=" << Info.TruncationError()
                << std::endl;
   }

   TruncErrSum += Info.TruncationError();

   // Update the effective Hamiltonian.
   HamR.push_front(contract_from_right(herm(*HPrev), *CPrev, HamR.front(), herm(*CPrev)));

   // Evolve the current site backwards in time.
   Iter = MaxIter;
   Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void
TDVP::EvolveLeftmostSite2(std::complex<double> Tau)
{
   // Form two-site centre block.
   LinearWavefunction::iterator CPrev = C;
   --CPrev;
   StateComponent C2 = local_tensor_prod(*CPrev, *C);

   BasicTriangularMPO::const_iterator HPrev = H;
   --HPrev;
   OperatorComponent H2 = local_tensor_prod(*HPrev, *H);

   HamL.pop_back();

   // Evolve two-site centre block.
   int Iter = MaxIter;
   double Err = ErrTol;

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, -I*Tau, Err, LogAmplitude);

   // Perform SVD on new C2.
   AMatSVD SL(C2, Tensor::ProductBasis<BasisList, BasisList>(CPrev->LocalBasis(), C->LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationErrorAbsolute(SL.begin(), SL.end(), SInfo, Info);
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, *CPrev, D, *C);
   *C = prod(D, *C);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Sites=" << Site-1 << "," << Site
                << " C2 Iter=" << Iter
                << " Err=" << Err
                << " States=" << Info.KeptStates()
		<< " TruncErr=" << Info.TruncationError()
                << std::endl;
   }

   MaxStates = Info.KeptStates();
   TruncErrSum += Info.TruncationError();

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*HPrev, herm(*CPrev), HamL.back(), *CPrev));
}

void
TDVP::IterateRight2(std::complex<double> Tau)
{
   // Evolve the current site backwards in time.
   int Iter = MaxIter;
   double Err = ErrTol;

   *C = LanczosExponential(*C, HEff1(HamL.back(), *H, HamR.front()), Iter, I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " C Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }

   // Form two-site centre block.
   ++Site;

   LinearWavefunction::iterator CPrev = C;
   ++C;
   StateComponent C2 = local_tensor_prod(*CPrev, *C);

   BasicTriangularMPO::const_iterator HPrev = H;
   ++H;
   OperatorComponent H2 = local_tensor_prod(*HPrev, *H);

   HamR.pop_front();

   // Evolve two-site centre block.
   Iter = MaxIter;
   Err = ErrTol;

   C2 = LanczosExponential(C2, HEff1(HamL.back(), H2, HamR.front()), Iter, -I*Tau, Err, LogAmplitude);

   // Perform SVD on new C2.
   AMatSVD SL(C2, Tensor::ProductBasis<BasisList, BasisList>(CPrev->LocalBasis(), C->LocalBasis()));
   TruncationInfo Info;
   AMatSVD::const_iterator Cutoff = TruncateFixTruncationErrorAbsolute(SL.begin(), SL.end(), SInfo, Info);
   RealDiagonalOperator D;
   SL.ConstructMatrices(SL.begin(), Cutoff, *CPrev, D, *C);
   *C = prod(D, *C);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Sites=" << Site-1 << "," << Site
                << " C2 Iter=" << Iter
                << " Err=" << Err
                << " States=" << Info.KeptStates()
		<< " TruncErr=" << Info.TruncationError()
                << std::endl;
   }

   MaxStates = std::max(MaxStates, Info.KeptStates());
   TruncErrSum += Info.TruncationError();

   // Update the effective Hamiltonian.
   HamL.push_back(contract_from_left(*HPrev, herm(*CPrev), HamL.back(), *CPrev));
}

void
TDVP::SweepLeft2(std::complex<double> Tau)
{
   while (Site > LeftStop + 1)
      this->IterateLeft2(Tau);

   this->EvolveLeftmostSite2(Tau);
}

void
TDVP::SweepRight2(std::complex<double> Tau)
{
   this->EvolveLeftmostSite2(Tau);

   while (Site < RightStop)
      this->IterateRight2(Tau);
}

void
TDVP::Evolve2()
{
   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   TruncErrSum = 0.0;

   std::vector<double>::const_iterator Alpha = Comp.Alpha.cbegin();
   std::vector<double>::const_iterator Beta = Comp.Beta.cbegin();

   while (Alpha != Comp.Alpha.cend())
   {
      this->UpdateHamiltonianLeft(Time, (*Alpha)*Timestep);

      this->SweepLeft2((*Alpha)*Timestep);
      Time += (*Alpha)*Timestep;
      ++Alpha;

      this->UpdateHamiltonianRight(Time, (*Beta)*Timestep);
      if (Ham.is_time_dependent())
      {
         ++H;
         HamR.pop_front();
      }

      this->SweepRight2((*Beta)*Timestep);
      Time += (*Beta)*Timestep;
      ++Beta;
   }

   if (Epsilon)
      this->CalculateEps();
}

void
TDVP::CalculateEps()
{
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   // Handle corner case where we have a one-site MPS (e.g. an IBC with a one-site window).
   if (LeftStop == RightStop)
      return;

   // Create a local copy so we can perform a single right-to-left sweep
   // without having to go back to the right.
   LinearWavefunction PsiLocal = Psi;
   LinearWavefunction::iterator CLocal = PsiLocal.end();
   // Move CLocal to the current position of C (may not be the final site).
   for (auto CCopy = C; CCopy != Psi.end(); ++CCopy)
      --CLocal;
   BasicTriangularMPO::const_iterator HLocal = H;
   std::deque<StateComponent> HamLLocal = HamL;
   std::deque<StateComponent> HamRLocal = HamR;
   int SiteLocal = RightStop;

   // Perform SVD to left-orthogonalize current site for NullSpace2.
   StateComponent CLeftOrtho = *CLocal;
   OrthogonalizeBasis2(CLeftOrtho);

   // Calculate error measure epsilon_1 and add to sum.
   StateComponent X = contract_from_left(*HLocal, herm(NullSpace2(CLeftOrtho)), HamLLocal.back(), *CLocal);
   double Eps1Sq = norm_frob_sq(scalar_prod(X, herm(HamRLocal.front())));
   Eps1SqSum += Eps1Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << SiteLocal
                << " Eps1Sq=" << Eps1Sq
                << std::endl;
   }

   while (SiteLocal > LeftStop)
   {
      // Perform SVD to right-orthogonalize current site.
      MatrixOperator U;
      RealDiagonalOperator D;
      std::tie(U, D) = OrthogonalizeBasis1(*CLocal);

      // Calculate the right half of epsilon_2 (see below).
      StateComponent Y = contract_from_right(herm(*HLocal), NullSpace1(*CLocal), HamRLocal.front(), herm(*CLocal));

      // Update the effective Hamiltonian.
      HamRLocal.push_front(contract_from_right(herm(*HLocal), *CLocal, HamRLocal.front(), herm(*CLocal)));

      // Move to the next site.
      --SiteLocal;
      --HLocal;
      --CLocal;

      *CLocal = prod(*CLocal, U*D);

      HamLLocal.pop_back();

      // Perform SVD to left-orthogonalize current site for NullSpace2.
      CLeftOrtho = *CLocal;
      OrthogonalizeBasis2(CLeftOrtho);

      // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
      StateComponent X = contract_from_left(*HLocal, herm(NullSpace2(CLeftOrtho)), HamLLocal.back(), *CLocal);
      Eps1Sq = norm_frob_sq(scalar_prod(X, herm(HamRLocal.front())));
      double Eps2Sq = norm_frob_sq(scalar_prod(X, herm(Y)));
      Eps1SqSum += Eps1Sq;
      Eps2SqSum += Eps2Sq;

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << SiteLocal
                   << " Eps1Sq=" << Eps1Sq
                   << " Eps2Sq=" << Eps2Sq
                   << std::endl;
      }
   }
}

void
TDVP::UpdateHamiltonianLeft(std::complex<double> t, std::complex<double> dt)
{
   if (!Ham.is_time_dependent())
      return;

   HamMPO = Ham(t, dt);
   H = HamMPO.end();
   --H;

   if (Verbose > 1)
      std::cout << "Recalculating left Hamiltonian environment..." << std::endl;

   HamL = std::deque<StateComponent>();
   HamL.push_back(Initial_E(HamMPO, Psi.Basis1()));

   LinearWavefunction::iterator CLocal = Psi.begin();
   BasicTriangularMPO::const_iterator HLocal = HamMPO.begin();
   int SiteLocal = LeftStop;

   while (SiteLocal < RightStop)
   {
      if (Verbose > 1)
         std::cout << "Site " << SiteLocal << std::endl;
      HamL.push_back(contract_from_left(*HLocal, herm(*CLocal), HamL.back(), *CLocal));
      ++HLocal, ++CLocal, ++SiteLocal;
   }
}

void
TDVP::UpdateHamiltonianRight(std::complex<double> t, std::complex<double> dt)
{
   if (!Ham.is_time_dependent())
      return;

   HamMPO = Ham(t, dt);
   H = HamMPO.begin();

   if (Verbose > 1)
      std::cout << "Recalculating right Hamiltonian environment..." << std::endl;

   HamR = std::deque<StateComponent>();
   HamR.push_front(Initial_F(HamMPO, Psi.Basis2()));

   LinearWavefunction::iterator CLocal = Psi.end();
   BasicTriangularMPO::const_iterator HLocal = HamMPO.end();
   int SiteLocal = RightStop + 1;

   while (SiteLocal > LeftStop + 1)
   {
      --HLocal, --CLocal, --SiteLocal;
      if (Verbose > 1)
         std::cout << "Site " << SiteLocal << std::endl;
      HamR.push_front(contract_from_right(herm(*HLocal), *CLocal, HamR.front(), herm(*CLocal)));
   }
}
