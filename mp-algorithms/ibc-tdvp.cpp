// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
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
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "ibc-tdvp.h"
#include "triangular_mpo_solver.h"
#include "tensor/regularize.h"
#include "tensor/tensor_eigen.h"
#include "lattice/unitcell-parser.h"
#include "linearalgebra/eigen.h"

WindowHamiltonian::WindowHamiltonian(std::string HamStrBackground, std::string HamStrWindow,
                                     int Size_, std::string Magnus_, std::string TimeVar_, int Verbose_)
   : Hamiltonian(HamStrBackground, Size_, Magnus_, TimeVar_, Verbose_)
{
   if (HamStrWindow.empty())
   {
      WindowEmpty = true;
      WindowTimeDependent = false;
   }
   else
   {
      WindowEmpty = false;
      InfiniteLattice WindowLattice;
      // NOTE: We cannot handle time-dependent windows at the moment. (TODO)
      //std::tie(WindowOperator, WindowLattice) = ParseUniCellOperatorStringAndLattice(HamStrWindow);
      std::tie(WindowMPO, WindowLattice) = ParseUnitCellOperatorAndLattice(HamStrWindow);
      WindowTimeDependent = false;

      // Cannot do this since operator== is not defined for InfiniteLattice
      // objects: we just assume that this is satified.
      //CHECK_EQUAL(Lattice, LatticeWindow);
   }
}

BasicTriangularMPO
WindowHamiltonian::operator()(int Left, int Right, std::complex<double> t, std::complex<double> dt) const
{
   // Initialize the background Hamiltonian for the whole unit cells inside the window.
   BasicTriangularMPO BackgroundMPO = this->Hamiltonian::operator()(t, dt);

   int UCSize = BackgroundMPO.size();

   // Find the left and right indices inside the window which are integer
   // multiples of the unit cell size.
   int LeftUC = -numerics::divp(-Left, UCSize).quot;
   int RightUC = numerics::divp(Right, UCSize).quot;
   // The remainders.
   int LeftRem = numerics::divp(-Left, UCSize).rem;
   int RightRem = numerics::divp(Right, UCSize).rem;

   // This happens when left and right are both inside of a single unit cell.
   // TODO: Can we handle this?
   if (RightUC < LeftUC)
      throw std::runtime_error("WindowHamiltonian: The desired range has no full unit cell. You must expand the window size manually to contain the operator in full unit cells.");

   // Extend the background Hamiltonian to fill the desired range with respect
   // to full unit cells (the leftover terms for partial unit cells are added
   // later).
   BasicTriangularMPO Result = repeat(BackgroundMPO, RightUC-LeftUC);

   // Add the window Hamiltonian if it exists.
   if (!WindowEmpty)
   {
      // TODO: We should incorporate extra sites into the window beforehand if this is not satisfied.
      if ((LeftUC * UCSize > WindowMPO.offset())
       || (RightUC * UCSize < WindowMPO.offset() + WindowMPO.size()))
         throw std::runtime_error("WindowHamiltonian: The range needs to contain the support of the window operators in full unit cells. You must expand the window size manually.");

      // We have to put the window MPO in triangular form so it can be added
      // onto the triangular background MPO.
      Result = Result + make_triangular(ExtendToCover(WindowMPO, (RightUC-LeftUC) * UCSize, LeftUC * UCSize));
   }

   // Now add the leftover terms to Hamiltonian by partially incorporating unit
   // cells from the left/right boundaries.
   std::vector<OperatorComponent> HamiltonianNew(Result.data().begin(), Result.data().end());

   auto HLeft = BackgroundMPO.end();
   for (int i = 0; i < LeftRem; ++i)
   {
      --HLeft;
      HamiltonianNew.insert(HamiltonianNew.begin(), *HLeft);
   }

   auto HRight = BackgroundMPO.begin();
   for (int i = 0; i < RightRem; ++i)
   {
      HamiltonianNew.insert(HamiltonianNew.end(), *HRight);
      ++HRight;
   }

   Result = BasicTriangularMPO(HamiltonianNew);

   return Result;
}

IBC_TDVP::IBC_TDVP(IBCWavefunction const& Psi_, WindowHamiltonian const& Ham_, IBC_TDVPSettings const& Settings_)
   : TDVP(Ham_, Settings_),
   HamWindow(Ham_), GMRESTol(Settings_.GMRESTol), FidTol(Settings_.FidTol), NExpand(Settings_.NExpand),
   Comoving(Settings_.Comoving), PsiLeft(Psi_.left()), PsiRight(Psi_.right())
{
   // We do not (currently) support time dependent Hamiltonians.
   // (We would need to recalculate the left and right block Hamiltonians each timestep.)
   CHECK(Ham.is_time_dependent() == false);

   //-----------------------------------
   // Handle the boundary wavefunctions.
   //-----------------------------------

   LeftQShift = Psi_.left_qshift();
   RightQShift = Psi_.right_qshift();
   WindowLeftSites = Psi_.window_left_sites();
   WindowRightSites = Psi_.window_right_sites();

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << '\n';

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
      std::cout << "Left energy = " << LeftEnergy << '\n';

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
   std::complex<double> BondEnergy = inner_prod(prod(PsiLeft.lambda_l(), prod(BlockHamL, PsiLeft.lambda_l())),
                                                delta_shift(BlockHamLR, PsiLeft.qshift()));

   if (Verbose > 0)
      std::cout << "Bond energy = " << BondEnergy << '\n';

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
      std::cout << "Right energy = " << RightEnergy << '\n';

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

   //-------------------
   // Handle the window.
   //-------------------

   Offset = Psi_.window_offset();
   LeftStop = Settings_.EvolutionWindowLeft;
   RightStop = Settings_.EvolutionWindowRight;
   Site = Offset;

   // Check whether there are any sites in the window, and if not, add some.
   if (Psi_.window_size() == 0)
   {
      // TODO: Expand window automatically.
      if (!HamWindow.is_window_empty())
         throw std::runtime_error("fatal: If using a window operator, your initial wavefunction window cannot be empty. You must expand the window size manually to contain the operator in full unit cells.");

      Psi = LinearWavefunction();
      MatrixOperator Lambda = Psi_.window().lambda_r();
      Lambda = Psi_.window().LeftU() * Lambda * Psi_.window().RightU();

      HamiltonianWindow = BasicTriangularMPO();

      if (Verbose > 1)
         std::cout << "Initial window size = 0, adding site to window..." << '\n';

      this->ExpandWindowRight();
      if (LeftStop == RightStop + 1)
         ++RightStop;

      // Incorporate the Lambda matrix into the wavefunction.
      *C = prod(Lambda, *C);

      HamR.pop_front();
   }
   else // Handle windows with > 0 sites.
   {
      MatrixOperator Lambda;
      std::tie(Psi, Lambda) = get_left_canonical(Psi_.window());

      C = Psi.begin();

      // Get the Hamiltonian for the window, being the sum of the background
      // Hamiltonian and the window Hamiltonian.
      HamiltonianWindow = HamWindow(Psi_.window_offset(), Psi_.window_offset() + Psi_.window_size());

      // Update iterators for the boundary MPOs to be added into the window when expanding.
      for (int i = 0; i < WindowLeftSites; ++i)
         --HLeft;

      for (int i = 0; i < WindowRightSites; ++i)
         ++HRight;

      H = HamiltonianWindow.begin();

      // Construct the left Hamiltonian environments for the window.
      while (C != Psi.end())
      {
         if (Verbose > 1)
            std::cout << "Site " << Site << '\n';
         HamL.push_back(contract_from_left(*H, herm(*C), HamL.back(), *C));
         MaxStates = std::max(MaxStates, (*C).Basis2().total_dimension());
         ++H, ++C, ++Site;
      }

      --H, --C, --Site;
      *C = prod(*C, Lambda);
      HamL.pop_back();
   }

   // Set the initial LeftStop for a comoving window.
   if (Comoving != 0)
      LeftStop = RightStop - Comoving + 1;

   // Ensure that the window contains LeftStop and RightStop.
   while (LeftStop < Offset)
      this->ExpandWindowLeft();

   while (RightStop > Offset + Psi.size()-1)
      this->ExpandWindowRight();

   //----------------------------------------------------
   // Initialize the evolution window expansion criteria.
   //----------------------------------------------------

   // Right-orthogonalize the window.
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

   // Add an extra site to the left of the evolution window if there isn't one.
   if (LeftStop == Offset)
      this->ExpandWindowLeft();

   // Save the left reference A-matrices.
   CRefLeft = *C;
   --C;
   CRefLeft2 = *C;
   ++C;

   // Left-orthogonalize the window.
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

   // Add an extra site to the right of the evolution window if there isn't one.
   if (RightStop == Offset + Psi.size()-1)
      this->ExpandWindowRight();

   // Save the right reference A-matrices.
   CRefRight = *C;
   ++C;
   CRefRight2 = *C;
   --C;
}

IBCWavefunction
IBC_TDVP::Wavefunction() const
{
   MatrixOperator I = MatrixOperator::make_identity(Psi.Basis2());
   WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Psi), I, Verbose-1);

   return IBCWavefunction(PsiLeft, PsiWindow, PsiRight, LeftQShift, RightQShift, Offset, WindowLeftSites, WindowRightSites);
}

void
IBC_TDVP::ExpandWindowLeft()
{
   --HLeft, --CLeft, --HamLeft;

   // Add the site to the window.
   Psi.push_front(delta_shift(*CLeft, LeftQShift));

   --Offset;
   ++WindowLeftSites;

   // Add the site's MPO to the window MPO.
   std::vector<OperatorComponent> HamiltonianNew;
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());
   HamiltonianNew.insert(HamiltonianNew.begin(), *HLeft);
   HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

   // Add the site's E matrix to the deque.
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
IBC_TDVP::ExpandWindowRight()
{
   // Add the site to the window.
   Psi.push_back(delta_shift(*CRight, RightQShift));

   ++WindowRightSites;

   // Add the site's MPO to the window MPO.
   std::vector<OperatorComponent> HamiltonianNew;
   HamiltonianNew.insert(HamiltonianNew.begin(), *HRight);
   HamiltonianNew.insert(HamiltonianNew.begin(), HamiltonianWindow.data().begin(), HamiltonianWindow.data().end());
   HamiltonianWindow = BasicTriangularMPO(HamiltonianNew);

   // Add the site's F matrix to the deque.
   HamR.push_back(*HamRight);

   // Reset iterators to previous location.
   C = Psi.end();
   H = HamiltonianWindow.end();
   for (int i = Offset + Psi.size()-1; i >= Site; --i)
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
IBC_TDVP::ExpandEvolutionWindowLeft()
{
   PRECONDITION(Site == LeftStop);

   // Expand the window if needed.
   while (LeftStop <= Offset + 1)
      this->ExpandWindowLeft();

   --LeftStop;

   // Obtain the left reference matrices.
   // We need to incorporate Lambda into CRefLeft, so we get that first.
   StateComponent CRightOrtho = *C;
   MatrixOperator U;
   RealDiagonalOperator D;

   std::tie(U, D) = OrthogonalizeBasis1(CRightOrtho);

   --C;
   CRefLeft = (*C)* U*D;
   --C;
   CRefLeft2 = *C;
   ++C;
   ++C;
}

void
IBC_TDVP::ExpandEvolutionWindowRight()
{
   PRECONDITION(Site == RightStop);

   // Expand the window if needed.
   while (RightStop >= Offset + Psi.size()-1 - 1)
      this->ExpandWindowRight();

   ++RightStop;

   // For a comoving window, update LeftStop as well.
   // NOTE: This will break CRefLeft(2), but we won't need them anyway,
   // so it should be fine(?)
   if (Comoving != 0)
      ++LeftStop;

   // Obtain the right reference matrices.
   // We need to incorporate Lambda into CRefRight, so we get that first.
   StateComponent CLeftOrtho = *C;
   MatrixOperator Vh;
   RealDiagonalOperator D;

   std::tie(D, Vh) = OrthogonalizeBasis2(CLeftOrtho);

   ++C;
   CRefRight = D*Vh * (*C);
   ++C;
   CRefRight2 = *C;
   --C;
   --C;
}

double
IBC_TDVP::CalculateFidelityLossLeft()
{
   --C;
   MatrixOperator Result = scalar_prod(herm(CRefLeft2), *C);
   ++C;
   Result = operator_prod(herm(CRefLeft), Result, *C);

   // The fidelity is given by the sum of singular values, so we calculate the SVD.
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SingularValueDecomposition(Result, U, D, Vh);

   return 1.0 - trace(D);
}

double
IBC_TDVP::CalculateFidelityLossRight()
{
   ++C;
   MatrixOperator Result = scalar_prod(*C, herm(CRefRight2));
   --C;
   Result = operator_prod(*C, Result, herm(CRefRight));

   // The fidelity is given by the sum of singular values, so we calculate the SVD.
   MatrixOperator U, Vh;
   RealDiagonalOperator D;
   SingularValueDecomposition(Result, U, D, Vh);

   return 1.0 - trace(D);
}

void
IBC_TDVP::SweepLeftEW(std::complex<double> Tau, bool Expand)
{
   while (Site > LeftStop)
   {
      if (Expand)
         this->ExpandLeft();
      this->EvolveCurrentSite(Tau);
      this->IterateLeft(Tau);
   }

   if (Expand)
   {
      // Add an extra site to the window if there isn't one.
      // (This shouldn't be able to happen anyway(?), so this is here just in case.)
      if (LeftStop == Offset)
         this->ExpandWindowLeft();
      this->ExpandLeft();
   }
   this->EvolveCurrentSite(Tau);

   // Expand the evolution window until the fidelity loss is below tolerance.
   double FidLoss = this->CalculateFidelityLossLeft();
   while (FidLoss > FidTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossLeft=" << FidLoss
                   << ", expanding window..." << '\n';

      this->ExpandEvolutionWindowLeft();

      this->IterateLeft(Tau);

      while (Site > LeftStop)
      {
         if (Expand)
            this->ExpandLeft();
         this->EvolveCurrentSite(Tau);
         this->IterateLeft(Tau);
      }

      if (Expand)
      {
         if (LeftStop == Offset)
            this->ExpandWindowLeft();
         this->ExpandLeft();
      }
      this->EvolveCurrentSite(Tau);

      FidLoss = this->CalculateFidelityLossLeft();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossLeft=" << FidLoss << '\n';
}

void
IBC_TDVP::SweepRightEW(std::complex<double> Tau, bool Expand)
{
   while (Site < RightStop)
   {
      if (Expand)
         this->ExpandRight();
      this->EvolveCurrentSite(Tau);
      this->IterateRight(Tau);
   }

   if (Expand)
   {
      // Add an extra site to the window if there isn't one.
      // (This shouldn't be able to happen anyway(?), so this is here just in case.)
      if (RightStop == Offset + Psi.size()-1)
         this->ExpandWindowRight();
      this->ExpandRight();
   }
   this->EvolveCurrentSite(Tau);

   // Expand the evolution window until the fidelity loss is below tolerance.
   double FidLoss = this->CalculateFidelityLossRight();
   while (FidLoss > FidTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossRight=" << FidLoss
                   << ", expanding window..." << '\n';

      this->ExpandEvolutionWindowRight();

      this->IterateRight(Tau);

      while (Site < RightStop)
      {
         if (Expand)
            this->ExpandRight();
         this->EvolveCurrentSite(Tau);
         this->IterateRight(Tau);
      }

      if (Expand)
      {
         if (RightStop == Offset + Psi.size()-1)
            this->ExpandWindowRight();
         this->ExpandRight();
      }
      this->EvolveCurrentSite(Tau);

      FidLoss = this->CalculateFidelityLossRight();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossRight=" << FidLoss << '\n';
}

void
IBC_TDVP::SweepRightFinalEW(std::complex<double> Tau, bool Expand)
{
   // Add an extra site to the window if there isn't one.
   // (This shouldn't be able to happen anyway(?), so this is here just in case.)
   if (RightStop == Offset + Psi.size()-1)
      this->ExpandWindowRight();

   if (Expand)
      this->ExpandRight();
   this->EvolveCurrentSite(Tau);
   this->CalculateEps1();

   while (Site < RightStop)
   {
      this->IterateRight(Tau);
      if (Expand)
         this->ExpandRight();
      this->EvolveCurrentSite(Tau);
      this->CalculateEps12();
   }

   double FidLoss = this->CalculateFidelityLossRight();
   while (FidLoss > FidTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossRight=" << FidLoss
                   << ", expanding window..." << '\n';

      this->ExpandEvolutionWindowRight();

      this->IterateRight(Tau);

      while (Site < RightStop)
      {
         if (Expand)
            this->ExpandRight();
         this->EvolveCurrentSite(Tau);
         this->CalculateEps12();
         this->IterateRight(Tau);
      }

      if (Expand)
      {
         if (RightStop == Offset + Psi.size()-1)
            this->ExpandWindowRight();
         this->ExpandRight();
      }
      this->EvolveCurrentSite(Tau);
      this->CalculateEps12();

      FidLoss = this->CalculateFidelityLossRight();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossRight=" << FidLoss << '\n';
}

void
IBC_TDVP::Evolve(bool Expand)
{
   ++TStep;
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Alpha = Comp.Alpha.cbegin();
   std::vector<double>::const_iterator Beta = Comp.Beta.cbegin();

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandEvolutionWindowRight();

   if (Comoving == 0)
      this->SweepLeftEW((*Alpha)*Timestep, Expand);
   else
      this->SweepLeft((*Alpha)*Timestep, Expand);
   ++Alpha;

   if (NExpand != 0 && Comoving == 0)
      if (TStep % NExpand == 0)
         this->ExpandEvolutionWindowLeft();

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
   {
      if (TStep % NExpand == 0)
      {
         // Expand the window if needed.
         while (RightStop >= Offset + Psi.size()-1 - 1)
            this->ExpandWindowRight();

         ++RightStop;

         // For a comoving window, update LeftStop as well.
         if (Comoving != 0)
            ++LeftStop;
      }
   }


   this->SweepLeft2((*Alpha)*Timestep);
   ++Alpha;

   if (NExpand != 0 && Comoving == 0)
   {
      if (TStep % NExpand == 0)
      {
         // Expand the window if needed.
         while (LeftStop <= Offset + 1)
            this->ExpandWindowLeft();

         --LeftStop;
      }
   }

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
