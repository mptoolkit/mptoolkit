// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-tdvp.cpp
//
// Copyright (C) 2021-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "ibc-tdvp.h"
#include "triangular_mpo_solver.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "linearalgebra/eigen.h"

IBC_TDVP::IBC_TDVP(IBCWavefunction const& Psi_, Hamiltonian const& Ham_, IBC_TDVPSettings const& Settings_)
   : TDVP(Ham_, Settings_),
   GMRESTol(Settings_.GMRESTol), FidTol(Settings_.FidTol), LambdaTol(Settings_.LambdaTol),
   UCExpand(Settings_.UCExpand), NExpand(Settings_.NExpand), Comoving(Settings_.Comoving),
   PsiLeft(Psi_.left()), PsiRight(Psi_.right())
{
   // We do not (currently) support time dependent Hamiltonians.
   // (We would need to recalculate the left and right block Hamiltonians each timestep.)
   CHECK(Ham.is_time_dependent() == false);

   LeftQShift = Psi_.left_qshift();
   RightQShift = Psi_.right_qshift();
   WindowLeftSites = Psi_.window_left_sites();
   WindowRightSites = Psi_.window_right_sites();

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   // Set left/right Hamiltonian sizes to match the unit cell sizes.
   HamiltonianLeft = Ham();
   HamiltonianRight = Ham();

   if (HamiltonianLeft.size() < PsiLeft.size())
      HamiltonianLeft = repeat(HamiltonianLeft, PsiLeft.size() / HamiltonianLeft.size());
   CHECK_EQUAL(HamiltonianLeft.size(), PsiLeft.size());

   if (HamiltonianRight.size() < PsiRight.size())
      HamiltonianRight = repeat(HamiltonianRight, PsiRight.size() / HamiltonianRight.size());
   CHECK_EQUAL(HamiltonianRight.size(), PsiRight.size());

   // Construct left Hamiltonian environment.
   StateComponent BlockHamL = Initial_E(HamiltonianLeft, PsiLeft.Basis1());
   std::complex<double> LeftEnergy = SolveHamiltonianMPO_Left(BlockHamL, PsiLeft, HamiltonianLeft,
                                                              GMRESTol, Verbose-1);

   if (Verbose > 0)
      std::cout << "Left energy = " << LeftEnergy << std::endl;

   // Remove a spurious contribution from the "bond energy", which is the
   // energy contribution from the terms in the Hamiltonian which cross the
   // bond at a unit cell boundary. To calculate this, we need the right
   // Hamiltonian environment for PsiLeft.
   LinearWavefunction PsiLinear;
   RealDiagonalOperator D;
   std::tie(D, PsiLinear) = get_right_canonical(PsiLeft);

   StateComponent BlockHamLR = Initial_F(HamiltonianLeft, PsiLinear.Basis2());
   MatrixOperator Rho = scalar_prod(D, herm(D));
   Rho = delta_shift(Rho, adjoint(PsiLeft.qshift()));

   SolveHamiltonianMPO_Right(BlockHamLR, PsiLinear, PsiLeft.qshift(), HamiltonianLeft, Rho, GMRESTol, Verbose-1);
   std::complex<double> BondEnergy = inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamLR);

   if (Verbose > 0)
      std::cout << "Bond energy = " << BondEnergy << std::endl;

   BlockHamL.back() -= BondEnergy * BlockHamL.front();

   // After calculating the bond energy, we apply the left qshift to the E matrix.
   BlockHamL = delta_shift(BlockHamL, LeftQShift);

   // Calculate the left Hamiltonian environments for each position in the unit cell.
   HamLeftUC.push_back(BlockHamL);

   HLeft = HamiltonianLeft.begin();
   CLeft = PsiLeft.begin();
   while (CLeft != PsiLeft.end())
   {
      StateComponent CShift = delta_shift(*CLeft, LeftQShift);
      HamLeftUC.push_back(contract_from_left(*HLeft, herm(CShift), HamLeftUC.back(), CShift));
      MaxStates = std::max(MaxStates, CShift.Basis2().total_dimension());
      ++HLeft, ++CLeft;
   }

   // Construct right Hamiltonian environment.
   StateComponent BlockHamR = Initial_F(HamiltonianRight, PsiRight.Basis2());
   std::complex<double> RightEnergy = SolveHamiltonianMPO_Right(BlockHamR, PsiRight, HamiltonianRight,
                                                                GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Right energy = " << RightEnergy << std::endl;

   BlockHamR = delta_shift(BlockHamR, RightQShift);

   // Calculate the right Hamiltonian environments for each position in the unit cell.
   HamRightUC.push_front(BlockHamR);

   HRight = HamiltonianRight.end();
   CRight = PsiRight.end();
   while (CRight != PsiRight.begin())
   {
      --HRight, --CRight;
      StateComponent CShift = delta_shift(*CRight, RightQShift);
      HamRightUC.push_front(contract_from_right(herm(*HRight), CShift, HamRightUC.front(), herm(CShift)));
      MaxStates = std::max(MaxStates, CShift.Basis1().total_dimension());
   }

   // Initialize the boundaries for the window Hamiltonian environment.
   Offset = Psi_.window_offset();
   RightStop = Psi_.window_size() - 1 + Psi_.window_offset();
   if (Comoving == 0)
      LeftStop = Offset;
   else
      LeftStop = RightStop - Comoving + 1;
   Site = Offset;

   HamLeft = HamLeftUC.cend();
   --HamLeft;
   for (int i = 0; i < WindowLeftSites; ++i)
      --CLeft, --HamLeft;

   // Note that we could just use *HamLeft if WindowLeftSites == 0, but this
   // would add a contribution to the energy from the boundary unit cell.
   if (WindowLeftSites == 0)
      HamL.push_back(delta_shift(BlockHamL, adjoint(PsiLeft.qshift())));
   else
      HamL.push_back(*HamLeft);

   HamRight = HamRightUC.cbegin();
   for (int i = 0; i < WindowRightSites; ++i)
      ++CRight, ++HamRight;

   if (WindowRightSites == 0)
      HamR.push_front(delta_shift(BlockHamR, PsiRight.qshift()));
   else
      HamR.push_front(*HamRight);

   ++HamRight;

   // Check whether there are any sites in the window, and if not, add some.
   if (Psi_.window_size() == 0)
   {
      Psi = LinearWavefunction();
      MatrixOperator Lambda = Psi_.window().lambda_r();
      Lambda = Psi_.window().LeftU() * Lambda * Psi_.window().RightU();

      HamiltonianWindow = BasicTriangularMPO();

      if (UCExpand)
      {
         if (Verbose > 1)
            std::cout << "Initial window size = 0, adding two unit cells to window..." << std::endl;

         this->ExpandWindowLeft();

         this->ExpandWindowRight();

         // Incorporate the Lambda matrix into the wavefunction.
         *C = prod(Lambda, *C);

         HamR.pop_front();
      }
      else
      {
         if (Verbose > 1)
            std::cout << "Initial window size = 0, adding site to window..." << std::endl;

         this->ExpandWindowRight();

         // Incorporate the Lambda matrix into the wavefunction.
         *C = prod(Lambda, *C);

         HamR.pop_front();
      }
   }
   else
   {
      MatrixOperator Lambda;
      std::tie(Psi, Lambda) = get_left_canonical(Psi_.window());

      C = Psi.begin();

      // The window size minus the sites partially incorporated from the
      // left/right boundaries.
      int WindowSizeWholeUCs = Psi.size() - WindowRightSites - WindowLeftSites;

      HamiltonianWindow = repeat(Ham(), WindowSizeWholeUCs / Ham().size());

      // Now add the "left over" terms to Hamiltonian from partially
      // incorporating unit cells from the left/right boundaries.
      std::vector<OperatorComponent> HamiltonianNew;
      HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());

      for (int i = 0; i < WindowLeftSites; ++i)
      {
         --HLeft;
         HamiltonianNew.insert(HamiltonianNew.begin(), *HLeft);
      }

      for (int i = 0; i < WindowRightSites; ++i)
      {
         HamiltonianNew.insert(HamiltonianNew.end(), *HRight);
         ++HRight;
      }

      HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

      H = HamiltonianWindow.begin();

      // Construct the left Hamiltonian environments for the window.
      while (C != Psi.end())
      {
         if (Verbose > 1)
            std::cout << "Site " << Site << std::endl;
         HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));
         MaxStates = std::max(MaxStates, (*C).Basis2().total_dimension());
         ++H, ++C, ++Site;
      }

      --H, --C, --Site;
      *C = prod(*C, Lambda);
      HamL.pop_back();

      if (UCExpand)
      {
         if (Verbose > 1 && (WindowLeftSites != 0 || WindowRightSites != 0))
            std::cout << "Expanding window to ensure boundary unit cell is fully incorporated into window..." << std::endl;

         while (WindowLeftSites != 0)
            this->ExpandWindowLeftSite();
         while (WindowRightSites != 0)
            this->ExpandWindowRightSite();
      }
   }

   // For a comoving window, ensure that the initial window has the correct
   // number of sites.
   if (Comoving != 0)
      while (LeftStop < Offset)
         this->ExpandWindowLeft();

   // Left-orthogonalize the window (should only be necessary if
   // UCExpand == true and we added extra sites on the right).
   while (Site < RightStop)
   {
      // Perform SVD to left-orthogonalize current site.
      MatrixOperator Vh;
      RealDiagonalOperator D;

      std::tie(D, Vh) = OrthogonalizeBasis2(*C);

      // Update the effective Hamiltonian.
      HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));

      // Move to the next site.
      ++Site;
      ++H;
      ++C;

      *C = prod(D*Vh, *C);

      HamR.pop_front();
   }
}

IBCWavefunction
IBC_TDVP::Wavefunction() const
{
   MatrixOperator I = MatrixOperator::make_identity(Psi.Basis2());
   WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Psi), I, Verbose-1);

   return IBCWavefunction(PsiLeft, PsiWindow, PsiRight, LeftQShift, RightQShift, Offset, WindowLeftSites, WindowRightSites);
}

void
IBC_TDVP::ExpandWindowLeftUC()
{
   // Add the unit cell to the window.
   auto CAdd = PsiLeft.end();
   while (CAdd != PsiLeft.begin())
   {
      --CAdd;
      Psi.push_front(delta_shift(*CAdd, LeftQShift));
   }

   // Change the leftmost index.
   Offset -= PsiLeft.size();
   if (Comoving == 0)
      LeftStop -= PsiLeft.size();

   // Add the unit cell to the Hamiltonian.
   std::vector<OperatorComponent> HamiltonianNew;
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianLeft.data().begin(), HamiltonianLeft.data().end());
   HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

   // Add the unit cell to the Hamiltonian environment.
   HamL.pop_front();
   HamL.insert(HamL.begin(), HamLeftUC.begin(), HamLeftUC.end());

   // Reset iterators to previous location.
   C = Psi.begin();
   H = HamiltonianWindow.begin();
   for (int i = Offset; i < Site; ++i)
      ++C, ++H;

   // Shift the quantum number of the boundary unit cell.
   LeftQShift = delta_shift(LeftQShift, PsiLeft.qshift());

   for (StateComponent& I : HamLeftUC)
      I = delta_shift(I, PsiLeft.qshift());
}

void
IBC_TDVP::ExpandWindowRightUC()
{
   // Add the unit cell to the window.
   auto CAdd = PsiRight.begin();
   while (CAdd != PsiRight.end())
   {
      Psi.push_back(delta_shift(*CAdd, RightQShift));
      ++CAdd;
   }

   // Change the rightmost index.
   RightStop += PsiRight.size();
   if (Comoving != 0)
      LeftStop += PsiRight.size();

   // Add the unit cell to the Hamiltonian.
   std::vector<OperatorComponent> HamiltonianNew;
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianRight.data().begin(), HamiltonianRight.data().end());
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());
   HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

   // Add the unit cell to the Hamiltonian environment.
   HamR.pop_back();
   HamR.insert(HamR.end(), HamRightUC.begin(), HamRightUC.end());

   // Reset iterators to previous location.
   C = Psi.end();
   H = HamiltonianWindow.end();
   for (int i = RightStop; i >= Site; --i)
      --C, --H;

   // Shift the quantum number of the boundary unit cell.
   RightQShift = delta_shift(RightQShift, adjoint(PsiRight.qshift()));

   for (StateComponent& I : HamRightUC)
      I = delta_shift(I, adjoint(PsiRight.qshift()));
}

void
IBC_TDVP::ExpandWindowLeftSite()
{
   --HLeft, --CLeft, --HamLeft;

   // Add the site to the window.
   Psi.push_front(delta_shift(*CLeft, LeftQShift));

   --Offset;
   if (Comoving == 0)
      --LeftStop;
   ++WindowLeftSites;

   // Add the site's operator to the Hamiltonian.
   std::vector<OperatorComponent> HamiltonianNew;
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());
   HamiltonianNew.insert(HamiltonianNew.begin(), *HLeft);
   HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

   // Add the site to the Hamiltonian environment.
   HamL.push_front(*HamLeft);

   // Reset iterators to previous location.
   C = Psi.begin();
   H = HamiltonianWindow.begin();
   for (int i = Offset; i < Site; ++i)
      ++C, ++H;

   if (CLeft == PsiLeft.begin())
   {
      // Shift the quantum number of the boundary unit cell.
      LeftQShift = delta_shift(LeftQShift, PsiLeft.qshift());

      for (StateComponent& I : HamLeftUC)
         I = delta_shift(I, PsiLeft.qshift());

      HLeft = HamiltonianLeft.end();
      CLeft = PsiLeft.end();
      HamLeft = HamLeftUC.cend();
      --HamLeft;

      WindowLeftSites = 0;
   }
}

void
IBC_TDVP::ExpandWindowRightSite()
{
   // Add the site to the window.
   Psi.push_back(delta_shift(*CRight, RightQShift));

   ++RightStop;
   if (Comoving != 0)
      ++LeftStop;
   ++WindowRightSites;

   // Add the site's operator to the Hamiltonian.
   std::vector<OperatorComponent> HamiltonianNew;
   HamiltonianNew.insert(HamiltonianNew.begin(), *HRight);
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());
   HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

   // Add the site to the Hamiltonian environment.
   HamR.push_back(*HamRight);

   // Reset iterators to previous location.
   C = Psi.end();
   H = HamiltonianWindow.end();
   for (int i = RightStop; i >= Site; --i)
      --C, --H;

   ++HRight, ++CRight, ++HamRight;
   if (CRight == PsiRight.end())
   {
      // Shift the quantum number of the boundary unit cell.
      RightQShift = delta_shift(RightQShift, adjoint(PsiRight.qshift()));

      for (StateComponent& I : HamRightUC)
         I = delta_shift(I, adjoint(PsiRight.qshift()));

      HRight = HamiltonianRight.begin();
      CRight = PsiRight.begin();
      HamRight = HamRightUC.cbegin();
      ++HamRight;

      WindowRightSites = 0;
   }
}


void
IBC_TDVP::ExpandWindowLeft()
{
   if (UCExpand)
      this->ExpandWindowLeftUC();
   else
      this->ExpandWindowLeftSite();
}
void
IBC_TDVP::ExpandWindowRight()
{
   if (UCExpand)
      this->ExpandWindowRightUC();
   else
      this->ExpandWindowRightSite();
}

double
IBC_TDVP::CalculateFidelityLossLeft()
{
   StateComponent CL;
   if (CLeft == PsiLeft.end())
      CL = delta_shift(PsiLeft[0], adjoint(PsiLeft.qshift()));
   else
      CL = *CLeft;
   CL = delta_shift(CL, LeftQShift);

   return 1.0 - norm_frob_sq(scalar_prod(herm(CL), *C));
}

double
IBC_TDVP::CalculateFidelityLossRight()
{
   StateComponent CR;
   if (CRight == PsiRight.begin())
      CR = delta_shift(PsiRight.get_back(), PsiRight.qshift());
   else
   {
      --CRight;
      CR = *CRight;
      ++CRight;
   }
   CR = delta_shift(CR, RightQShift);

   return 1.0 - norm_frob_sq(scalar_prod(*C, herm(CR)));
}

double
IBC_TDVP::CalculateLambdaDiffLeft()
{
   // Right-orthogonalize current site to find LambdaL.
   StateComponent CRightOrtho = *C;
   MatrixOperator U;
   RealDiagonalOperator D;

   std::tie(U, D) = OrthogonalizeBasis1(CRightOrtho);

   // Ensure that the left and right bases of LambdaL are the same.
   MatrixOperator LambdaL = (U*D)*herm(U);

   MatrixOperator LambdaLeft = PsiLeft.lambda(PsiLeft.size() - WindowLeftSites);
   LambdaLeft = delta_shift(LambdaLeft, LeftQShift);

   return norm_frob_sq(LambdaL - LambdaLeft);
}

double
IBC_TDVP::CalculateLambdaDiffRight()
{
   // Left-orthogonalize current site to find LambdaR.
   StateComponent CLeftOrtho = *C;
   MatrixOperator Vh;
   RealDiagonalOperator D;

   std::tie(D, Vh) = OrthogonalizeBasis2(CLeftOrtho);

   // Ensure that the left and right bases of LambdaR are the same.
   MatrixOperator LambdaR = herm(Vh)*(D*Vh);

   MatrixOperator LambdaRight = PsiRight.lambda(WindowRightSites);
   LambdaRight = delta_shift(LambdaRight, RightQShift);

   return norm_frob_sq(LambdaR - LambdaRight);
}

// TODO: At the moment, the handling of expanding the bond at the edge of the
// window is inefficient: we first evolve the leftmost site and then check if
// we want to expand the window, and if we do, we unevolve the site, expand the
// bond, then evolve again.
// It would be better to try to expand the bond before evolving the leftmost
// site, and if we do need to expand the bond, we can just expand the window as
// well.
void
IBC_TDVP::SweepLeftEW(std::complex<double> Tau, bool Expand)
{
   while (Site > LeftStop)
   {
      if (Expand && C->Basis1().total_dimension() < SInfo.MaxStates)
         this->ExpandLeft();
      this->EvolveCurrentSite(Tau);
      this->IterateLeft(Tau);
   }

   StateComponent Tmp = *C;
   this->EvolveCurrentSite(Tau);

   double FidLoss = this->CalculateFidelityLossLeft();
   double LambdaDiff = this->CalculateLambdaDiffLeft();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossLeft=" << FidLoss
                   << " LambdaDiffLeft=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowLeft();

      if (Expand && C->Basis1().total_dimension() < SInfo.MaxStates)
      {
         *C = Tmp;
         this->ExpandLeft();
         this->EvolveCurrentSite(Tau);
      }
      this->IterateLeft(Tau);

      while (Site > LeftStop)
      {
         if (Expand && C->Basis1().total_dimension() < SInfo.MaxStates)
            this->ExpandLeft();
         this->EvolveCurrentSite(Tau);
         this->IterateLeft(Tau);
      }

      Tmp = *C;
      this->EvolveCurrentSite(Tau);

      FidLoss = this->CalculateFidelityLossLeft();
      LambdaDiff = this->CalculateLambdaDiffLeft();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossLeft=" << FidLoss
                << " LambdaDiffLeft=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::SweepRightEW(std::complex<double> Tau, bool Expand)
{
   while (Site < RightStop)
   {
      if (Expand && C->Basis2().total_dimension() < SInfo.MaxStates)
         this->ExpandRight();
      this->EvolveCurrentSite(Tau);
      this->IterateRight(Tau);
   }

   StateComponent Tmp = *C;
   this->EvolveCurrentSite(Tau);

   double FidLoss = this->CalculateFidelityLossRight();
   double LambdaDiff = this->CalculateLambdaDiffRight();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossRight=" << FidLoss
                   << " LambdaDiffRight=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowRight();

      if (Expand && C->Basis2().total_dimension() < SInfo.MaxStates)
      {
         *C = Tmp;
         this->ExpandRight();
         this->EvolveCurrentSite(Tau);
      }
      this->IterateRight(Tau);

      while (Site < RightStop)
      {
         if (Expand && C->Basis2().total_dimension() < SInfo.MaxStates)
            this->ExpandRight();
         this->EvolveCurrentSite(Tau);
         this->IterateRight(Tau);
      }

      Tmp = *C;
      this->EvolveCurrentSite(Tau);

      FidLoss = this->CalculateFidelityLossRight();
      LambdaDiff = this->CalculateLambdaDiffRight();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossRight=" << FidLoss
                << " LambdaDiffRight=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::SweepRightFinalEW(std::complex<double> Tau, bool Expand)
{
   StateComponent Tmp = *C;

   if (Expand && Site < RightStop && C->Basis2().total_dimension() < SInfo.MaxStates)
      this->ExpandRight();
   this->EvolveCurrentSite(Tau);
   this->CalculateEps1();

   while (Site < RightStop)
   {
      this->IterateRight(Tau);
      Tmp = *C;
      if (Expand && Site < RightStop && C->Basis2().total_dimension() < SInfo.MaxStates)
         this->ExpandRight();
      this->EvolveCurrentSite(Tau);
      this->CalculateEps12();
   }

   double FidLoss = this->CalculateFidelityLossRight();
   double LambdaDiff = this->CalculateLambdaDiffRight();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossRight=" << FidLoss
                   << " LambdaDiffRight=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowRight();

      if (Expand && C->Basis2().total_dimension() < SInfo.MaxStates)
      {
         *C = Tmp;
         this->ExpandRight();
         this->EvolveCurrentSite(Tau);
      }
      this->IterateRight(Tau);

      while (Site < RightStop)
      {
         if (Expand && C->Basis2().total_dimension() < SInfo.MaxStates)
            this->ExpandRight();
         this->EvolveCurrentSite(Tau);
         this->CalculateEps12();
         this->IterateRight(Tau);
      }

      Tmp = *C;
      this->EvolveCurrentSite(Tau);
      this->CalculateEps12();

      FidLoss = this->CalculateFidelityLossRight();
      LambdaDiff = this->CalculateLambdaDiffRight();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossRight=" << FidLoss
                << " LambdaDiffRight=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::Evolve(bool Expand)
{
   ++TStep;
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Alpha = Comp.Alpha.cbegin();
   std::vector<double>::const_iterator Beta = Comp.Beta.cbegin();

   if (NExpand != 0 && Comoving == 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowLeft();

   if (Comoving == 0)
      this->SweepLeftEW((*Alpha)*Timestep, Expand);
   else
      this->SweepLeft((*Alpha)*Timestep, Expand);
   ++Alpha;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowRight();

   while (Alpha != Comp.Alpha.cend())
   {
      this->SweepRightEW((*Beta)*Timestep, Expand);
      ++Beta;

      if (Comoving == 0)
         this->SweepLeftEW((*Alpha)*Timestep, Expand);
      else
         this->SweepLeft((*Alpha)*Timestep, Expand);
      ++Alpha;
   }

   if (Epsilon)
      this->SweepRightFinalEW((*Beta)*Timestep, Expand);
   else
      this->SweepRightEW((*Beta)*Timestep, Expand);
}

void
IBC_TDVP::Evolve2()
{
   ++TStep;
   TruncErrSum = 0.0;

   std::vector<double>::const_iterator Alpha = Comp.Alpha.cbegin();
   std::vector<double>::const_iterator Beta = Comp.Beta.cbegin();

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowLeft();

   this->SweepLeft2((*Alpha)*Timestep);
   ++Alpha;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowRight();

   this->SweepRight2((*Beta)*Timestep);
   ++Beta;

   while (Alpha != Comp.Alpha.cend())
   {
      this->SweepLeft2((*Alpha)*Timestep);
      ++Alpha;

      this->SweepRight2((*Beta)*Timestep);
      ++Beta;
   }

   if (Epsilon)
      this->CalculateEps();
}
