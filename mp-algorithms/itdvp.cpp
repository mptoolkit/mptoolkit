// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/itdvp.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021-2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
   std::complex<double> dt = Comp.Gamma.back()*Timestep;
   HamMPO = Ham(Time-dt, dt);

   // Make sure that Psi and HamMPO have the same unit cell.
   InfiniteWavefunctionLeft PsiCanonical = Psi_;
   int UnitCellSize = statistics::lcm(PsiCanonical.size(), HamMPO.size());

   if (PsiCanonical.size() != UnitCellSize)
   {
      std::cout << "Warning: Extending wavefunction unit cell to " << UnitCellSize << " sites." << std::endl;
      PsiCanonical = repeat(PsiCanonical, UnitCellSize / PsiCanonical.size());
      Ham.set_size(UnitCellSize);
      std::complex<double> dt = Comp.Gamma.back()*Timestep;
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
   MatrixOperator U;
   RealDiagonalOperator D;
   std::tie(U, D, PsiR) = get_right_canonical(PsiCanonical);

   PsiR.set_back(prod(PsiR.get_back(), delta_shift(U, adjoint(QShift))));

   BlockHamR = Initial_F(HamMPO, PsiR.Basis2());
   Rho = scalar_prod(D, herm(D));
   Rho = delta_shift(Rho, adjoint(QShift));
   SolveHamiltonianMPO_Right(BlockHamR, PsiR, QShift, HamMPO, Rho, GMRESTol, Verbose-1);

   U = delta_shift(U, adjoint(QShift));
   BlockHamR = prod(prod(U, BlockHamR), herm(U));

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
   MatrixOperator M = ExpandBasis1(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(U*Vh, *C);
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
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(*C, U*Vh);
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
   std::deque<StateComponent> HamLOld = HamL;
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
   std::deque<StateComponent> HamROld = HamR;
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
iTDVP::Evolve()
{
   Time = InitialTime + ((double) TStep)*Timestep;
   ++TStep;

   XYCalculated = false;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();

   while (Gamma != Comp.Gamma.cend())
   {
      // If we have already updated the Hamiltonian before expanding the bonds,
      // then we should not do it again.
      if (!HamUpdated)
         this->UpdateHamiltonianLeft(Time, (*Gamma)*Timestep);
      else
         HamUpdated = false;

      this->EvolveLeft((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;

      this->UpdateHamiltonianRight(Time, (*Gamma)*Timestep);

      this->EvolveRight((*Gamma)*Timestep);
      Time += (*Gamma)*Timestep;
      ++Gamma;
   }

   if (Epsilon)
      this->CalculateEps();
}

void
iTDVP::CalculateEps()
{
   X = std::deque<StateComponent>();
   Y = std::deque<StateComponent>();

   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   *C = prod(*C, LambdaR);

   // Right-orthogonalize the unit cell.
   while (Site > LeftStop)
   {
      // Perform SVD to right-orthogonalize current site.
      MatrixOperator M = ExpandBasis1(*C);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *C = prod(Vh, *C);

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
      MatrixOperator M = ExpandBasis1(*C);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      StateComponent CRightOrtho = prod(Vh, *C);
      *C = prod(U*D*Vh, *C);

      // Calculate the right half of epsilon_2 for the left end of the unit cell.
      Y.push_back(contract_from_right(herm(*H), NullSpace1(CRightOrtho), HamR.front(), herm(*C)));
   }

   while (Site < RightStop)
   {
      // Perform SVD to left-orthogonalize current site.
      MatrixOperator M = ExpandBasis2(*C);
      MatrixOperator U, Vh;
      RealDiagonalOperator D;

      SingularValueDecomposition(M, U, D, Vh);

      *C = prod(*C, U);

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

      if (Epsilon)
      {
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
   }

   E = inner_prod(contract_from_left(*H, herm(*C), HamL.back(), *C), HamR.front());

   // Perform SVD to left-orthogonalize the rightmost site.
   MatrixOperator M = ExpandBasis2(*C);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   *C = prod(*C, U*Vh);
   LambdaR = herm(Vh)*(D*Vh);

   // Update left block Hamiltonian.
   BlockHamL = contract_from_left(*H, herm(*C), HamL.back(), *C);
   BlockHamL.back() -= E * BlockHamL.front();
   HamL.front() = delta_shift(BlockHamL, QShift);

   // Calculate the left half of epsilon_2.
   X.push_back(contract_from_left(*H, herm(NullSpace2(*C)), HamL.back(), *C));

   Y.push_back(delta_shift(Y.front(), adjoint(QShift)));
   Y.pop_front();

   if (Epsilon)
   {
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
         this->CalculateEpsN();
   }

   XYCalculated = true;
}

void
iTDVP::CalculateEpsN()
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

void
iTDVP::ExpandRightBond()
{
   if ((*C).Basis2().total_dimension() < SInfo.MaxStates)
   {
      // Take the truncated SVD of P_2 H|Psi>.
      CMatSVD SL(scalar_prod(X.front(), herm(Y.front())));
      TruncationInfo Info;
      StatesInfo SInfoLocal = SInfo;
      // Subtract the current bond dimension from the number of additional states to be added.
      SInfoLocal.MinStates = std::max(0, SInfoLocal.MinStates - (*C).Basis2().total_dimension());
      SInfoLocal.MaxStates = std::max(0, SInfoLocal.MaxStates - (*C).Basis2().total_dimension());

      CMatSVD::const_iterator Cutoff;
      if (ForceExpand)
      {
         Cutoff = SL.begin();
         for (int i = 0; Cutoff != SL.end() && i < SInfoLocal.MaxStates; ++i)
            ++Cutoff;
      }
      else
         Cutoff = TruncateFixTruncationError(SL.begin(), SL.end(), SInfoLocal, Info);

      MatrixOperator U, Vh;
      RealDiagonalOperator D;
      SL.ConstructMatrices(SL.begin(), Cutoff, U, D, Vh);

      // Construct new basis.
      SumBasis<VectorBasis> NewBasis((*C).Basis2(), U.Basis2());
      // Construct a unitary to regularize the new basis.
      MatrixOperator UReg = Regularize(NewBasis);

      MaxStates = std::max(MaxStates, NewBasis.total_dimension());

      // Add the new states to CExpand.
      StateComponent CExpand = tensor_row_sum(*COld, prod(NullSpace2(*COld), U), NewBasis);
      CExpand = prod(CExpand, herm(UReg));

      // Add a zero tensor of the same dimensions to C.
      StateComponent Z = StateComponent((*C).LocalBasis(), (*C).Basis1(), U.Basis2());
      *C = tensor_row_sum(*C, Z, NewBasis);
      *C = prod(*C, herm(UReg));

      if (Verbose > 1)
      {
         std::cout << "Timestep=" << TStep
                   << " Site=" << Site
                   << " NewStates=" << Info.KeptStates()
                   << " TotalStates=" << NewBasis.total_dimension()
                   << std::endl;
      }

      // Move to next site, handling the rightmost site separately.
      if (Site < RightStop)
      {
         HamL.push_back(contract_from_left(*H, herm(CExpand), HamLOld.front(), CExpand));

         HamLOld.pop_front();
         ++C;
         ++COld;
      }
      else
      {
         BlockHamL = contract_from_left(*H, herm(CExpand), HamLOld.front(), CExpand);
         BlockHamL.back() -= E * BlockHamL.front();
         HamL.front() = delta_shift(BlockHamL, QShift);

         C = Psi.begin();
         *C = delta_shift(*C, adjoint(QShift));

         // Add zeros to LambdaR so the bonds match.
         MatrixOperator Z = MatrixOperator(Vh.Basis1(), LambdaR.Basis2());
         LambdaR = tensor_col_sum(LambdaR, Z, NewBasis);
         LambdaR = UReg * LambdaR;
      }

      ++Site;
      ++H;
      X.pop_front();
      Y.pop_front();

      // Add zeros to the current site so the left bond matches.
      Z = StateComponent((*C).LocalBasis(), Vh.Basis1(), (*C).Basis2());
      *C = tensor_col_sum(*C, Z, NewBasis);
      *C = prod(UReg, *C);
   }
   else // If we don't need to add any more states.
   {
      // Move to next site, handling the rightmost site separately.
      if (Site < RightStop)
      {
         HamL.push_back(contract_from_left(*H, herm(*COld), HamLOld.front(), *COld));

         HamLOld.pop_front();
         ++C;
         ++COld;
      }
      else
      {
         BlockHamL = contract_from_left(*H, herm(*COld), HamLOld.front(), *COld);
         BlockHamL.back() -= E * BlockHamL.front();
         HamL.front() = delta_shift(BlockHamL, QShift);

         C = Psi.begin();
         *C = delta_shift(*C, adjoint(QShift));
      }

      ++Site;
      ++H;
      X.pop_front();
      Y.pop_front();
   }
}

void
iTDVP::ExpandBonds()
{
   if (!HamUpdated)
      this->UpdateHamiltonianLeft(Time, Comp.Gamma.front()*Timestep);
   HamUpdated = true;

   if (!XYCalculated)
      this->CalculateEps();

   C = Psi.begin();
   H = HamMPO.begin();
   Site = LeftStop;
   HamLOld = HamL;
   HamL = std::deque<StateComponent>(1, BlockHamL);

   LinearWavefunction PsiOld = Psi;
   COld = PsiOld.begin();

   while (Site <= RightStop)
      this->ExpandRightBond();

   *C = delta_shift(*C, QShift);

   C = Psi.end();
   --C;
   H = HamMPO.end();
   --H;
   Site = RightStop;

   XYCalculated = false;
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

   MatrixOperator Rho = scalar_prod(LambdaR, herm(LambdaR));
   Rho = delta_shift(Rho, QShift);

   BlockHamL = Initial_E(HamMPO, Psi.Basis1());
   BlockHamL.back() = HamL.front().back();

   SolveHamiltonianMPO_Left(BlockHamL, Psi, QShift, HamMPO, Rho, GMRESTol, Verbose-1);
   HamL = std::deque<StateComponent>(1, BlockHamL);
   BlockHamL = delta_shift(BlockHamL, adjoint(QShift));

   MatrixOperator RhoOld = scalar_prod(LambdaROld, herm(LambdaROld));
   RhoOld = delta_shift(RhoOld, adjoint(QShift));

   BlockHamR = Initial_F(HamMPO, PsiOld.Basis2());
   BlockHamR.front() = HamR.back().front();

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

   X = std::deque<StateComponent>();
   Y = std::deque<StateComponent>();
   XYCalculated = false;
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

   X = std::deque<StateComponent>();
   Y = std::deque<StateComponent>();
   XYCalculated = false;
}
