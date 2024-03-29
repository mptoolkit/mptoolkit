// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/itdvp.cpp
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

#include "itdvp.h"
#include "lanczos-exponential-new.h"
#include "triangular_mpo_solver.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"
#include "common/statistics.h"

std::complex<double> const I(0.0, 1.0);

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

iTDVP::iTDVP(InfiniteWavefunctionLeft const& Psi_, Hamiltonian const& Ham_, iTDVPSettings Settings_)
   : TDVP(Ham_, Settings_),
     GMRESTol(Settings_.GMRESTol), MaxSweeps(Settings_.MaxSweeps),
     LambdaTol(Settings_.LambdaTol), NEps(Settings_.NEps)
{
   // Initialize Psi and Ham.
   Time = InitialTime;
   std::complex<double> dt = Comp.Beta.back()*Timestep;
   HamMPO = Ham(Time-dt, dt);

   // Make sure that Psi and HamMPO have the same unit cell.
   InfiniteWavefunctionLeft PsiCanonical = Psi_;
   int UnitCellSize = statistics::lcm(PsiCanonical.size(), HamMPO.size());

   if (PsiCanonical.size() != UnitCellSize)
   {
      std::cout << "Warning: Extending wavefunction unit cell to " << UnitCellSize << " sites." << std::endl;
      PsiCanonical = repeat(PsiCanonical, UnitCellSize / PsiCanonical.size());
      Ham.set_size(UnitCellSize);
      std::complex<double> dt = Comp.Beta.back()*Timestep;
      HamMPO = Ham(Time-dt, dt);
   }

   QShift = PsiCanonical.qshift();
   LogAmplitude = PsiCanonical.log_amplitude();

   std::tie(Psi, LambdaR) = get_left_canonical(PsiCanonical);

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   H = HamMPO.begin();

   BlockHamL = Initial_E(HamMPO, Psi.Basis1());
   MatrixOperator Rho = scalar_prod(LambdaR, herm(LambdaR));
   Rho = delta_shift(Rho, QShift);
   InitialE = SolveHamiltonianMPO_Left(BlockHamL, Psi, QShift, HamMPO, Rho, GMRESTol, Verbose-1);
   HamL.push_back(BlockHamL);
   BlockHamL = delta_shift(BlockHamL, adjoint(QShift));

   for (auto const& I : Psi)
   {
      if (Verbose > 1)
         std::cout << "Site " << (HamL.size()) << std::endl;
      HamL.push_back(contract_from_left(*H, herm(I), HamL.back(), I));
      MaxStates = std::max(MaxStates, (I).Basis2().total_dimension());
      ++H;
   }

   // Calculate initial right Hamiltonian.
   LinearWavefunction PsiR;
   RealDiagonalOperator D;
   std::tie(D, PsiR) = get_right_canonical(PsiCanonical);

   BlockHamR = Initial_F(HamMPO, PsiR.Basis2());
   Rho = scalar_prod(D, herm(D));
   Rho = delta_shift(Rho, adjoint(QShift));
   SolveHamiltonianMPO_Right(BlockHamR, PsiR, QShift, HamMPO, Rho, GMRESTol, Verbose-1);

   HamR.push_front(BlockHamR);
   BlockHamR = delta_shift(BlockHamR, QShift);

   // Initialize to the right-most site.
   HamL.pop_back();
   C = Psi.end();
   --C;
   --H;
   Site = Psi.size() - 1;

   LeftStop = 0;
   RightStop = Psi.size() - 1;
}

InfiniteWavefunctionLeft
iTDVP::Wavefunction() const
{
   return InfiniteWavefunctionLeft::Construct(Psi, QShift, Normalize ? 0.0 : LogAmplitude);
}

void
iTDVP::OrthogonalizeLeftmostSite()
{
   E = inner_prod(HamL.back(), contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

   // Right-orthogonalize current site.
   MatrixOperator U;
   RealDiagonalOperator D;

   std::tie(U, D) = OrthogonalizeBasis1(*C);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(U, *C);
   LambdaR = (U*D)*herm(U);
   LambdaR = delta_shift(LambdaR, adjoint(QShift));

   // Update right block Hamiltonian.
   BlockHamR = contract_from_right(herm(*H), *C, HamR.front(), herm(*C));
   BlockHamR.front() -= E * BlockHamR.back();
   HamR.back() = delta_shift(BlockHamR, adjoint(QShift));
}

void
iTDVP::OrthogonalizeRightmostSite()
{
   E = inner_prod(contract_from_left(*H, herm(*C), HamL.back(), *C), HamR.front());

   // Left-orthogonalize current site.
   MatrixOperator Vh;
   RealDiagonalOperator D;

   std::tie(D, Vh) = OrthogonalizeBasis2(*C);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(*C, Vh);
   LambdaR = herm(Vh)*(D*Vh);
   LambdaR = delta_shift(LambdaR, QShift);

   // Update left block Hamiltonian.
   BlockHamL = contract_from_left(*H, herm(*C), HamL.back(), *C);
   BlockHamL.back() -= E * BlockHamL.front();
   HamL.front() = delta_shift(BlockHamL, QShift);
}

void
iTDVP::EvolveLambdaRRight(std::complex<double> Tau)
{
   int Iter = MaxIter;
   double Err = ErrTol;

   LambdaR = LanczosExponential(LambdaR, HEff2(BlockHamL, HamR.back()), Iter, I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " LambdaR Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void
iTDVP::EvolveLambdaRLeft(std::complex<double> Tau)
{
   int Iter = MaxIter;
   double Err = ErrTol;

   LambdaR = LanczosExponential(LambdaR, HEff2(HamL.front(), BlockHamR), Iter, I*Tau, Err, LogAmplitude);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " LambdaR Iter=" << Iter
                << " Err=" << Err
                << std::endl;
   }
}

void
iTDVP::EvolveLeft(std::complex<double> Tau)
{
   int Sweep = 0;
   double FidelityLoss = 1.0;
   double LambdaRDiff = 1.0;

   PsiOld = Psi;
   LambdaROld = LambdaR;
   HamLOld = HamL;
   double LogAmplitudeOld = LogAmplitude;

   do {
      ++Sweep;

      LinearWavefunction PsiPrev = Psi;
      MatrixOperator LambdaRPrev = LambdaR;

      Psi = PsiOld;
      C = Psi.end();
      --C;
      H = HamMPO.end();
      --H;
      Site = RightStop;

      HamL = HamLOld;
      HamR = std::deque<StateComponent>(1, delta_shift(BlockHamR, adjoint(QShift)));

      if (Sweep > 1)
         this->EvolveLambdaRRight(Tau);

      LogAmplitude = LogAmplitudeOld;

      *C = prod(*C, LambdaR);

      this->SweepLeft(Tau);

      this->OrthogonalizeLeftmostSite();

      if (Sweep > 1)
      {
         // Calculate the fidelity loss compared to the previous sweep.
         MatrixOperator Rho = scalar_prod(herm(LambdaRPrev), LambdaR);
         Rho = delta_shift(Rho, QShift);
         Rho = inject_left(Rho, Psi, PsiPrev);

         MatrixOperator U, Vh;
         RealDiagonalOperator D;
         SingularValueDecomposition(Rho, U, D, Vh);

         FidelityLoss = 1.0 - trace(D);

         LambdaRDiff = norm_frob_sq(LambdaR - LambdaRPrev);

         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepL=" << Sweep
                      << " FidelityLoss=" << FidelityLoss
                      << " LambdaRDiff=" << LambdaRDiff
                      << std::endl;
         }
      }
      else
      {
         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepL=" << Sweep
                      << std::endl;
         }
      }
   }
   while (LambdaRDiff > LambdaTol && Sweep < MaxSweeps);

   if (Sweep == MaxSweeps)
      std::cout << "WARNING: MaxSweeps reached, LambdaRDiff=" << LambdaRDiff << std::endl;

   LambdaR = delta_shift(LambdaR, QShift);
}

void
iTDVP::EvolveRight(std::complex<double> Tau)
{
   int Sweep = 0;
   double FidelityLoss = 1.0;
   double LambdaRDiff = 1.0;

   PsiOld = Psi;
   LambdaROld = LambdaR;
   HamROld = HamR;
   double LogAmplitudeOld = LogAmplitude;

   do {
      ++Sweep;

      LinearWavefunction PsiPrev = Psi;
      MatrixOperator LambdaRPrev = LambdaR;

      Psi = PsiOld;
      C = Psi.begin();
      H = HamMPO.begin();
      Site = LeftStop;

      HamL = std::deque<StateComponent>(1, delta_shift(BlockHamL, QShift));
      HamR = HamROld;

      if (Sweep > 1)
         this->EvolveLambdaRLeft(Tau);

      LogAmplitude = LogAmplitudeOld;

      *C = prod(LambdaR, *C);

      this->SweepRight(Tau);

      this->OrthogonalizeRightmostSite();

      if (Sweep > 1)
      {
         // Calculate the fidelity loss compared to the previous sweep.
         MatrixOperator Rho = scalar_prod(LambdaR, herm(LambdaRPrev));
         Rho = delta_shift(Rho, adjoint(QShift));
         Rho = inject_right(Rho, Psi, PsiPrev);

         MatrixOperator U, Vh;
         RealDiagonalOperator D;
         SingularValueDecomposition(Rho, U, D, Vh);

         FidelityLoss = 1.0 - trace(D);

         LambdaRDiff = norm_frob_sq(LambdaR - LambdaRPrev);

         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepR=" << Sweep
                      << " FidelityLoss=" << FidelityLoss
                      << " LambdaRDiff=" << LambdaRDiff
                      << std::endl;
         }
      }
      else
      {
         if (Verbose > 0)
         {
            std::cout << "Timestep=" << TStep
                      << " SweepR=" << Sweep
                      << std::endl;
         }
      }
   }
   while (LambdaRDiff > LambdaTol && Sweep < MaxSweeps);

   if (Sweep == MaxSweeps)
      std::cout << "WARNING: MaxSweeps reached, LambdaRDiff=" << LambdaRDiff << std::endl;

   LambdaR = delta_shift(LambdaR, adjoint(QShift));
}

void
iTDVP::Evolve(bool Expand)
{
   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   std::vector<double>::const_iterator Alpha = Comp.Alpha.cbegin();
   std::vector<double>::const_iterator Beta = Comp.Beta.cbegin();

   while (Alpha != Comp.Alpha.cend())
   {
      // We do not need/cannot update the Hamiltonian on the first timestep.
      if (TStep > 1)
         this->UpdateHamiltonianLeft(Time, (*Alpha)*Timestep);

      if (Expand)
         this->ExpandBondsLeft();

      this->EvolveLeft((*Alpha)*Timestep);
      Time += (*Alpha)*Timestep;
      ++Alpha;

      this->UpdateHamiltonianRight(Time, (*Beta)*Timestep);

      if (Expand)
         this->ExpandBondsRight();

      this->EvolveRight((*Beta)*Timestep);
      Time += (*Beta)*Timestep;
      ++Beta;
   }

   if (Epsilon)
      this->CalculateEps();
}

void
iTDVP::CalculateEps()
{
   std::deque<StateComponent> X, Y;

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   *C = prod(*C, LambdaR);

   // Right-orthogonalize the unit cell.
   while (Site > LeftStop)
   {
      // Perform SVD to right-orthogonalize current site.
      MatrixOperator U;
      RealDiagonalOperator D;

      std::tie(U, D) = OrthogonalizeBasis1(*C);

      // Update the effective Hamiltonian.
      HamR.push_front(contract_from_right(herm(*H), *C, HamR.front(), herm(*C)));

      // Move to the next site.
      --Site;
      --H;
      --C;

      *C = prod(*C, U*D);

      HamL.pop_back();
   }

   {
      // Perform SVD to right-orthogonalize current site for NullSpace1.
      StateComponent CRightOrtho = *C;
      OrthogonalizeBasis1(CRightOrtho);

      // Calculate the right half of epsilon_2 for the left end of the unit cell.
      Y.push_back(contract_from_right(herm(*H), NullSpace1(CRightOrtho), HamR.front(), herm(*C)));
   }

   while (Site < RightStop)
   {
      // Perform SVD to left-orthogonalize current site.
      MatrixOperator Vh;
      RealDiagonalOperator D;

      std::tie(D, Vh) = OrthogonalizeBasis2(*C);

      // Calculate the left half of epsilon_2.
      X.push_back(contract_from_left(*H, herm(NullSpace2(*C)), HamL.back(), *C));

      // Update the effective Hamiltonian.
      HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));

      // Move to the next site.
      ++Site;
      ++H;
      ++C;

      StateComponent CRightOrtho = prod(Vh, *C);
      *C = prod(D*Vh, *C);

      HamR.pop_front();

      // Calculate the right half of epsilon_2.
      Y.push_back(contract_from_right(herm(*H), NullSpace1(CRightOrtho), HamR.front(), herm(*C)));

      // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
      double Eps1Sq = norm_frob_sq(scalar_prod(HamL.back(), herm(Y.back())));
      double Eps2Sq = norm_frob_sq(scalar_prod(X.back(), herm(Y.back())));
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

   // We handle the rightmost site separately.
   E = inner_prod(contract_from_left(*H, herm(*C), HamL.back(), *C), HamR.front());

   // Perform SVD to left-orthogonalize the rightmost site.
   MatrixOperator Vh;
   RealDiagonalOperator D;

   std::tie(D, Vh) = OrthogonalizeBasis2(*C);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(*C, Vh);
   LambdaR = herm(Vh)*(D*Vh);

   // Update left block Hamiltonian.
   BlockHamL = contract_from_left(*H, herm(*C), HamL.back(), *C);
   BlockHamL.back() -= E * BlockHamL.front();
   HamL.front() = delta_shift(BlockHamL, QShift);

   // Calculate the left half of epsilon_2.
   X.push_back(contract_from_left(*H, herm(NullSpace2(*C)), HamL.back(), *C));

   Y.push_back(delta_shift(Y.front(), adjoint(QShift)));
   Y.pop_front();

   // Calculate error measures epsilon_1 and epsilon_2 and add to sums.
   double Eps1Sq = norm_frob_sq(scalar_prod(BlockHamL, herm(Y.back())));
   double Eps2Sq = norm_frob_sq(scalar_prod(X.back(), herm(Y.back())));
   Eps1SqSum += Eps1Sq;
   Eps2SqSum += Eps2Sq;

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << LeftStop
                << " Eps1Sq=" << Eps1Sq
                << " Eps2Sq=" << Eps2Sq
                << std::endl;
   }

   if (NEps > 2)
      this->CalculateEpsN(X, Y);
}

void
iTDVP::CalculateEpsN(std::deque<StateComponent> X, std::deque<StateComponent> Y)
{
   EpsNSqSum = std::vector<double>(NEps-2, 0.0);
   std::deque<StateComponent>::const_iterator Xi = X.end();
   --Xi;
   std::deque<StateComponent>::const_iterator Yi = Y.end();
   --Yi;

   while (Site >= LeftStop)
   {
      LinearWavefunction::iterator CiLocal = C;
      BasicTriangularMPO::const_iterator HLocal = H;
      std::deque<StateComponent>::const_iterator XiLocal = Xi;
      StateComponent YLocal = *Yi;
      int SiteLocal = Site;
      int NShifts = 0;

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site;
      }

      for (int i = 0; i < NEps-2; ++i)
      {
         StateComponent CLocal = *CiLocal;
         for (int j = 0; j < NShifts; ++j)
            CLocal = delta_shift(CLocal, QShift);

         StateComponent I = StateComponent::ConstructFullBasis1((CLocal).LocalBasis(), YLocal.Basis1());
         YLocal = contract_from_right(herm(*HLocal), I, YLocal, herm(CLocal));

         if (SiteLocal == LeftStop)
         {
            CiLocal = Psi.end();
            HLocal = HamMPO.end();
            XiLocal = X.end();
            SiteLocal = RightStop+1;
            ++NShifts;
         }
         --CiLocal, --HLocal, --XiLocal, --SiteLocal;

         StateComponent XLocal = *XiLocal;
         for (int j = 0; j < NShifts; ++j)
            XLocal = delta_shift(XLocal, QShift);

         double EpsNSq = norm_frob_sq(scalar_prod(XLocal, herm(YLocal)));
         EpsNSqSum[i] += EpsNSq;

         if (Verbose > 1)
            std::cout << " Eps" << i+3 << "Sq=" << EpsNSq;
      }

      if (Verbose > 1)
         std::cout << std::endl;

      --C, --H, --Xi, --Yi, --Site;
   }

   C = Psi.end();
   --C;
   H = HamMPO.end();
   --H;
   Site = RightStop;
}

// The aim here is to expand all of the bonds in the unit cell from right to
// left without modifying the actual wavefunction. But, for the expansion
// algorithm, we need the orthogonality center to be on right site of the bond
// being expanded. To solve this, we store the current orthogonality center in
// a class variable CCenter, while we keep the actual unit cell in left
// orthogonal form.
//
// Another issue arises from the unit cell boundary, since we want an
// unmodified version of the rightmost A-matrix when evolving the bond across
// the boundary. We save a copy of the rightmost A-matrix to deal with this.
void
iTDVP::ExpandBondsLeft()
{
   MaxStates = 0;

   HamLOld = HamL;
   HamLOld.push_front(delta_shift(HamLOld.back(), QShift));
   HamLOld.pop_back();
   HamL = std::deque<StateComponent>();

   HamROld = HamR;

   CBoundary = delta_shift(*C, QShift);
   CCenter = prod(*C, LambdaR);

   while (Site > LeftStop)
   {
      this->ExpandLeft();
      --C, --H, --Site;
   }

   // Handle leftmost site separately.
   this->ExpandLeft();

   // Move back to the right end of the unit cell.
   C = Psi.end();
   --C;
   H = HamMPO.end();
   --H;
   Site = RightStop;

   // Update the left block Hamiltonian.
   BlockHamL = delta_shift(HamL.front(), adjoint(QShift));
}

void
iTDVP::ExpandLeft()
{
   auto CNext = C;
   // Shift to the end of the unit cell if we are at the beginning.
   if (CNext == Psi.begin())
      CNext = Psi.end();
   --CNext;

   auto HNext = H;
   if (HNext == HamMPO.begin())
      HNext = HamMPO.end();
   --HNext;

   // We need to handle the last bond in the unit cell separately: see below.
   StateComponent CExpand = Site == LeftStop ? CBoundary : *CNext;

   VectorBasis AddBasis(C->Basis1().GetSymmetryList());
   int NewStates = 0;
   if (CExpand.Basis2().total_dimension() < SInfo.MaxStates)
   {
      TruncationInfo Info;
      std::tie(Info, AddBasis) = ExpandLeftEnvironment(CExpand, CCenter, HamLOld.back(), HamROld.front(), *HNext, *H, SInfo);
      NewStates = Info.KeptStates();
   }

   // We need to save this for dealing with the last bond.
   if (Site == RightStop)
      AddBasisBoundary = AddBasis;

   // Add the zeros to the left-orthogonal A matrix in the unit cell.
   SumBasis<VectorBasis> NewBasis(C->Basis1(), AddBasis);
   Regularizer R(NewBasis);
   StateComponent ZR(C->LocalBasis(), AddBasis, C->Basis2());
   *C = RegularizeBasis1(R, tensor_col_sum(*C, ZR, NewBasis));

   int TotalStates = CExpand.Basis2().total_dimension();
   MaxStates = std::max(MaxStates, TotalStates);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " NewStates=" << NewStates
                << " TotalStates=" << TotalStates
                << std::endl;
   }

   // Calculate updated E matrix with added states.
   HamL.push_front(contract_from_left(*HNext, herm(CExpand), HamLOld.back(), CExpand));
   HamLOld.pop_back();

   if (Site == LeftStop)
   {
      // We have to add the zeros to CExpand after doing the expansion for the
      // last bond, since adding the zeros in beforehand will result in a
      // larger null space tensor.
      CExpand = delta_shift(CExpand, adjoint(QShift));
      SumBasis<VectorBasis> NewBasisBoundary(CExpand.Basis1(), AddBasisBoundary);
      Regularizer RBoundary(NewBasisBoundary);
      StateComponent ZRBoundary(CExpand.LocalBasis(), AddBasisBoundary, CExpand.Basis2());
      *CNext = RegularizeBasis1(RBoundary, tensor_col_sum(CExpand, ZRBoundary, NewBasisBoundary));

      // Add zeros the lambda matrix so we can multiply it back onto the unit cell.
      LambdaR = delta_shift(LambdaR, QShift);
      MatrixOperator ZLambda = MatrixOperator(AddBasis, LambdaR.Basis2());
      LambdaR = RegularizeBasis1(R, tensor_col_sum(LambdaR, ZLambda, NewBasis));
      LambdaR = delta_shift(LambdaR, adjoint(QShift));
   }
   else
   {
      *CNext = CExpand;

      // Right orthogonalize current site.
      MatrixOperator U;
      RealDiagonalOperator D;
      std::tie(U, D) = OrthogonalizeBasis1(CCenter);

      // Calculate F matrix using new right-orthogonal A matrix.
      HamROld.push_front(contract_from_right(herm(*H), CCenter, HamROld.front(), herm(CCenter)));

      // Move the orthgonality center.
      CCenter = prod(CExpand, U*D);
   }
}

void
iTDVP::ExpandBondsRight()
{
   MaxStates = 0;

   HamROld = HamR;
   HamROld.push_back(delta_shift(HamROld.front(), adjoint(QShift)));
   HamROld.pop_front();
   HamR = std::deque<StateComponent>();

   HamLOld = HamL;

   CBoundary = delta_shift(*C, adjoint(QShift));
   CCenter = prod(LambdaR, *C);

   while (Site < RightStop)
   {
      this->ExpandRight();
      ++C, ++H, ++Site;
   }

   // Handle rightmost site separately.
   this->ExpandRight();

   // Move back to the left end of the unit cell.
   C = Psi.begin();
   H = HamMPO.begin();
   Site = LeftStop;

   // Update the right block Hamiltonian.
   BlockHamR = delta_shift(HamR.back(), QShift);
}

void
iTDVP::ExpandRight()
{
   auto CNext = C;
   ++CNext;
   // Shift to the beginning of the unit cell if we are at the end.
   if (CNext == Psi.end())
      CNext = Psi.begin();

   auto HNext = H;
   ++HNext;
   if (HNext == HamMPO.end())
      HNext = HamMPO.begin();

   // We need to handle the last bond in the unit cell separately: see below.
   StateComponent CExpand = Site == RightStop ? CBoundary : *CNext;

   VectorBasis AddBasis(C->Basis2().GetSymmetryList());
   int NewStates = 0;
   if (CExpand.Basis1().total_dimension() < SInfo.MaxStates)
   {
      TruncationInfo Info;
      std::tie(Info, AddBasis) = ExpandRightEnvironment(CCenter, CExpand, HamLOld.back(), HamROld.front(), *H, *HNext, SInfo);
      NewStates = Info.KeptStates();
   }

   // We need to save this for dealing with the last bond.
   if (Site == LeftStop)
      AddBasisBoundary = AddBasis;

   // Add the zeros to the right-orthogonal A matrix in the unit cell.
   SumBasis<VectorBasis> NewBasis(C->Basis2(), AddBasis);
   Regularizer R(NewBasis);
   StateComponent ZL(C->LocalBasis(), C->Basis1(), AddBasis);
   *C = RegularizeBasis2(tensor_row_sum(*C, ZL, NewBasis), R);

   int TotalStates = CExpand.Basis1().total_dimension();
   MaxStates = std::max(MaxStates, TotalStates);

   if (Verbose > 1)
   {
      std::cout << "Timestep=" << TStep
                << " Site=" << Site
                << " NewStates=" << NewStates
                << " TotalStates=" << TotalStates
                << std::endl;
   }

   // Calculate updated F matrix with added states.
   HamR.push_back(contract_from_right(herm(*HNext), CExpand, HamROld.front(), herm(CExpand)));
   HamROld.pop_front();

   if (Site == RightStop)
   {
      // We have to add the zeros to CExpand after doing the expansion for the
      // last bond, since adding the zeros in beforehand will result in a
      // larger null space tensor.
      CExpand = delta_shift(CExpand, QShift);
      SumBasis<VectorBasis> NewBasisBoundary(CExpand.Basis2(), AddBasisBoundary);
      Regularizer RBoundary(NewBasisBoundary);
      StateComponent ZLBoundary(CExpand.LocalBasis(), CExpand.Basis1(), AddBasisBoundary);
      *CNext = RegularizeBasis2(tensor_row_sum(CExpand, ZLBoundary, NewBasisBoundary), RBoundary);

      // Add zeros the lambda matrix so we can multiply it back onto the unit cell.
      LambdaR = delta_shift(LambdaR, adjoint(QShift));
      MatrixOperator ZLambda = MatrixOperator(LambdaR.Basis1(), AddBasis);
      LambdaR = RegularizeBasis2(tensor_row_sum(LambdaR, ZLambda, NewBasis), R);
      LambdaR = delta_shift(LambdaR, QShift);
   }
   else
   {
      *CNext = CExpand;

      // Left orthogonalize current site.
      MatrixOperator Vh;
      RealDiagonalOperator D;
      std::tie(D, Vh) = OrthogonalizeBasis2(CCenter);

      // Calculate E matrix using new left-orthogonal A matrix.
      HamLOld.push_back(contract_from_left(*H, herm(CCenter), HamLOld.back(), CCenter));

      // Move the orthgonality center.
      CCenter = prod(D*Vh, CExpand);
   }
}

void
iTDVP::UpdateHamiltonianLeft(std::complex<double> t, std::complex<double> dt)
{
   if (!Ham.is_time_dependent())
      return;

   HamMPO = Ham(t, dt);
   H = HamMPO.end();
   --H;

   if (Verbose > 1)
      std::cout << "Recalculating Hamiltonian environments (left)..." << std::endl;

   MatrixOperator Rho = scalar_prod(herm(LambdaR), LambdaR);
   Rho = delta_shift(Rho, QShift);

   BlockHamL = Initial_E(HamMPO, Psi.Basis1());
   BlockHamL.back() = HamL.front().back();

   SolveHamiltonianMPO_Left(BlockHamL, Psi, QShift, HamMPO, Rho, GMRESTol, Verbose-1);
   HamL = std::deque<StateComponent>(1, BlockHamL);
   BlockHamL = delta_shift(BlockHamL, adjoint(QShift));

   MatrixOperator RhoOld = scalar_prod(herm(LambdaROld), LambdaROld);
   RhoOld = delta_shift(RhoOld, adjoint(QShift));

   BlockHamR = Initial_F(HamMPO, PsiOld.Basis2());
   BlockHamR.front() = HamROld.back().front();

   SolveHamiltonianMPO_Right(BlockHamR, PsiOld, QShift, HamMPO, RhoOld, GMRESTol, Verbose-1);
   HamR = std::deque<StateComponent>(1, BlockHamR);
   BlockHamR = delta_shift(BlockHamR, QShift);

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
iTDVP::UpdateHamiltonianRight(std::complex<double> t, std::complex<double> dt)
{
   if (!Ham.is_time_dependent())
      return;

   HamMPO = Ham(t, dt);
   H = HamMPO.begin();

   if (Verbose > 1)
      std::cout << "Recalculating Hamiltonian environments (right)..." << std::endl;

   MatrixOperator RhoOld = scalar_prod(LambdaROld, herm(LambdaROld));
   RhoOld = delta_shift(RhoOld, QShift);

   BlockHamL = Initial_E(HamMPO, PsiOld.Basis1());
   BlockHamL.back() = HamL.front().back();

   SolveHamiltonianMPO_Left(BlockHamL, PsiOld, QShift, HamMPO, RhoOld, GMRESTol, Verbose-1);
   HamL = std::deque<StateComponent>(1, BlockHamL);
   BlockHamL = delta_shift(BlockHamL, adjoint(QShift));

   MatrixOperator Rho = scalar_prod(LambdaR, herm(LambdaR));
   Rho = delta_shift(Rho, adjoint(QShift));

   BlockHamR = Initial_F(HamMPO, Psi.Basis2());
   BlockHamR.front() = HamR.back().front();

   SolveHamiltonianMPO_Right(BlockHamR, Psi, QShift, HamMPO, Rho, GMRESTol, Verbose-1);
   HamR = std::deque<StateComponent>(1, BlockHamR);
   BlockHamR = delta_shift(BlockHamR, QShift);

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
