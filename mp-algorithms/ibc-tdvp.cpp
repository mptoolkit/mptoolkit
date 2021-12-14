// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-tdvp.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2021 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

IBC_TDVP::IBC_TDVP(IBCWavefunction const& Psi_, BasicTriangularMPO const& Ham_,
                   std::complex<double> Timestep_, Composition Comp_, int MaxIter_,
                   double ErrTol_, double GMRESTol_, double FidTol_, double LambdaTol_,
                   int NExpand_, StatesInfo SInfo_, int Verbose_)
   : TDVP(Ham_, Timestep_, Comp_, MaxIter_, ErrTol_, SInfo_, Verbose_),
   GMRESTol(GMRESTol_), FidTol(FidTol_), LambdaTol(LambdaTol_), NExpand(NExpand_)
{
   PsiLeft = Psi_.Left;
   PsiRight = Psi_.Right;

   if (Verbose > 0)
      std::cout << "Constructing Hamiltonian block operators..." << std::endl;

   // Set left/right Hamiltonian sizes to match the unit cell sizes.
   HamiltonianLeft = Hamiltonian;
   HamiltonianRight = Hamiltonian;

   if (HamiltonianLeft.size() < PsiLeft.size())
      HamiltonianLeft = repeat(HamiltonianLeft, PsiLeft.size() / HamiltonianLeft.size());

   if (HamiltonianRight.size() < PsiRight.size())
      HamiltonianRight = repeat(HamiltonianRight, PsiRight.size() / HamiltonianRight.size());

   // Construct left Hamiltonian environment.
   StateComponent BlockHamL = Initial_E(HamiltonianLeft, PsiLeft.Basis2());
   std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamiltonianLeft,
                                                         GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Starting energy (left eigenvalue) = " << LeftEnergy << std::endl;

   HamLeftL.push_back(BlockHamL);

   BasicTriangularMPO::const_iterator HLeft = HamiltonianLeft.begin();
   InfiniteWavefunctionLeft::const_mps_iterator CLeft = PsiLeft.begin();
   while (CLeft != PsiLeft.end())
   {
      HamLeftL.push_back(contract_from_left(*HLeft, herm(*CLeft), HamLeftL.back(), *CLeft));
      MaxStates = std::max(MaxStates, (*CLeft).Basis2().total_dimension());
      ++HLeft, ++CLeft;
   }

   // Construct right Hamiltonian environment.
   StateComponent BlockHamR = Initial_F(HamiltonianRight, PsiRight.Basis1());
   std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiRight, HamiltonianRight,
                                                           GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Starting energy (right eigenvalue) = " << RightEnergy << std::endl;

   HamRightR.push_front(BlockHamR);

   BasicTriangularMPO::const_iterator HRight = HamiltonianRight.end();
   InfiniteWavefunctionLeft::const_mps_iterator CRight = PsiRight.end();
   while (CRight != PsiRight.begin())
   {
      --HRight, --CRight;
      HamRightR.push_front(contract_from_right(herm(*HRight), *CRight, HamRightR.front(), herm(*CRight)));
      MaxStates = std::max(MaxStates, (*CRight).Basis1().total_dimension());
   }

   // Construct window Hamiltonian environment.
   LeftStop = Psi_.window_offset();
   RightStop = Psi_.window_size() - 1 + Psi_.window_offset();
   Site = LeftStop;

   HamL.push_back(BlockHamL);
   HamR.push_front(BlockHamR);

   // Check whether there are any sites in the window, and if not, add two unit cells.
   if (Psi_.window_size() == 0)
   {
      if (Verbose > 1)
         std::cout << "Initial window size = 0, adding two unit cells to window..." << std::endl;

      Psi = LinearWavefunction();
      MatrixOperator Lambda = Psi_.Window.lambda_r();
      Lambda = Psi_.Window.LeftU() * Lambda * Psi_.Window.RightU();

      Hamiltonian = BasicTriangularMPO();

      this->ExpandWindowLeft();

      this->ExpandWindowRight();

      // Incorporate the Lambda matrix into the wavefunction.
      *C = prod(Lambda, *C);

      HamR.pop_front();

      // Left-orthogonalize the window.
      while (Site < RightStop)
      {
         // Perform SVD to left-orthogonalize current site.
         MatrixOperator M = ExpandBasis2(*C);
         MatrixOperator U, Vh;
         RealDiagonalOperator D;

         SingularValueDecomposition(M, U, D, Vh);

         *C = prod(*C, U);

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
   else
   {
      MatrixOperator Lambda;
      std::tie(Psi, Lambda) = get_left_canonical(Psi_.Window);

      C = Psi.begin();

      if (Hamiltonian.size() < Psi.size())
         Hamiltonian = repeat(Hamiltonian, Psi.size() / Hamiltonian.size());
      H = Hamiltonian.begin();

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
   }
}

IBCWavefunction
IBC_TDVP::Wavefunction() const
{
   MatrixOperator I = MatrixOperator::make_identity(Psi.Basis2());
   WavefunctionSectionLeft PsiWindow = WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Psi), I, Verbose-1);

   return IBCWavefunction(PsiLeft, PsiWindow, PsiRight, LeftStop);
}

void
IBC_TDVP::ExpandWindowLeft()
{
   // Add the unit cell to the window.
   LinearWavefunction PsiCell;
   RealDiagonalOperator Lambda;
   std::tie(PsiCell, Lambda) = get_left_canonical(PsiLeft);

   Psi.push_front(PsiCell);

   // Add the unit cell to the Hamiltonian.
   std::vector<OperatorComponent> HamNew;
   HamNew.insert(HamNew.begin(), Hamiltonian.data().begin(), Hamiltonian.data().end());
   HamNew.insert(HamNew.begin(), HamiltonianLeft.data().begin(), HamiltonianLeft.data().end());
   Hamiltonian = BasicTriangularMPO(HamNew);

   // Add the unit cell to the Hamiltonian environment.
   HamL.pop_front();
   HamL.insert(HamL.begin(), HamLeftL.begin(), HamLeftL.end());

   // Change the leftmost index.
   LeftStop -= PsiLeft.size();

   // Reset iterators to previous location.
   C = Psi.begin();
   H = Hamiltonian.begin();
   for (int i = LeftStop; i < Site; ++i)
      ++C, ++H;
}

void
IBC_TDVP::ExpandWindowRight()
{
   // Add the unit cell to the window.
   LinearWavefunction PsiCell;
   RealDiagonalOperator Lambda;
   std::tie(Lambda, PsiCell) = get_right_canonical(PsiRight);

   Psi.push_back(PsiCell);

   // Add the unit cell to the Hamiltonian.
   std::vector<OperatorComponent> HamNew;
   HamNew.insert(HamNew.begin(), HamiltonianRight.data().begin(), HamiltonianRight.data().end());
   HamNew.insert(HamNew.begin(), Hamiltonian.data().begin(), Hamiltonian.data().end());
   Hamiltonian = BasicTriangularMPO(HamNew);

   // Add the unit cell to the Hamiltonian environment.
   HamR.pop_back();
   HamR.insert(HamR.end(), HamRightR.begin(), HamRightR.end());

   // Change the rightmost index.
   RightStop += PsiRight.size();

   // Reset iterators to previous location.
   C = Psi.begin();
   H = Hamiltonian.begin();
   for (int i = LeftStop; i < Site; ++i)
      ++C, ++H;
}

double
IBC_TDVP::CalculateFidelityLossLeft()
{
   StateComponent CL = PsiLeft[0];
   return 1.0 - norm_frob_sq(scalar_prod(herm(CL), *C));
}

double
IBC_TDVP::CalculateFidelityLossRight()
{
   StateComponent CR = PsiRight.get_back();
   return 1.0 - norm_frob_sq(scalar_prod(*C, herm(CR)));
}

double
IBC_TDVP::CalculateLambdaDiffLeft()
{
   // Right-orthogonalize current site to find LambdaL.
   StateComponent CCopy = *C;
   MatrixOperator M = ExpandBasis1(CCopy);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaL are the same.
   MatrixOperator LambdaL = (U*D)*herm(U);

   MatrixOperator LambdaLBoundary = PsiLeft.lambda_r();

   return norm_frob_sq(LambdaL - LambdaLBoundary);
}

double
IBC_TDVP::CalculateLambdaDiffRight()
{
   // Left-orthogonalize current site to find LambdaR.
   StateComponent CCopy = *C;
   MatrixOperator M = ExpandBasis2(CCopy);
   MatrixOperator U, Vh;
   RealDiagonalOperator D;

   SingularValueDecomposition(M, U, D, Vh);

   // Ensure that the left and right bases of LambdaR are the same.
   MatrixOperator LambdaR = herm(Vh)*(D*Vh);

   MatrixOperator LambdaRBoundary = PsiRight.lambda_l();

   return norm_frob_sq(LambdaR - LambdaRBoundary);
}

void
IBC_TDVP::SweepLeftEW(std::complex<double> Tau)
{
   this->SweepLeft(Tau);

   double FidLoss = this->CalculateFidelityLossLeft();
   double LambdaDiff = this->CalculateLambdaDiffLeft();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossLeft=" << FidLoss
                   << " LambdaDiffLeft=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowLeft();

      this->IterateLeft(Tau);
      this->SweepLeft(Tau);

      FidLoss = this->CalculateFidelityLossLeft();
      LambdaDiff = this->CalculateLambdaDiffLeft();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossLeft=" << FidLoss
                << " LambdaDiffLeft=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::SweepRightEW(std::complex<double> Tau)
{
   this->SweepRight(Tau);

   double FidLoss = this->CalculateFidelityLossRight();
   double LambdaDiff = this->CalculateLambdaDiffRight();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossRight=" << FidLoss
                   << " LambdaDiffRight=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowRight();

      this->IterateRight(Tau);
      this->SweepRight(Tau);

      FidLoss = this->CalculateFidelityLossRight();
      LambdaDiff = this->CalculateLambdaDiffRight();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossRight=" << FidLoss
                << " LambdaDiffRight=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::SweepRightFinalEW(std::complex<double> Tau)
{
   this->SweepRightFinal(Tau);

   double FidLoss = this->CalculateFidelityLossRight();
   double LambdaDiff = this->CalculateLambdaDiffRight();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossRight=" << FidLoss
                   << " LambdaDiffRight=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowRight();

      while (Site < RightStop)
      {
         this->IterateRight(Tau);
         this->EvolveCurrentSite(Tau);
         this->CalculateEps12();
      }
      FidLoss = this->CalculateFidelityLossRight();
      LambdaDiff = this->CalculateLambdaDiffRight();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossRight=" << FidLoss
                << " LambdaDiffRight=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::SweepLeftExpandEW(std::complex<double> Tau)
{
   this->SweepLeftExpand(Tau);

   double FidLoss = this->CalculateFidelityLossLeft();
   double LambdaDiff = this->CalculateLambdaDiffLeft();
   while (FidLoss > FidTol || LambdaDiff > LambdaTol)
   {
      if (Verbose > 0)
         std::cout << "FidelityLossLeft=" << FidLoss
                   << " LambdaDiffLeft=" << LambdaDiff
                   << ", expanding window..." << std::endl;

      this->ExpandWindowLeft();

      // TODO: Implement better handling of the expansion of bonds between the window and boundary.
      if ((*C).Basis1().total_dimension() < SInfo.MaxStates)
         this->ExpandLeftBond();
      this->IterateLeft(Tau);
      this->SweepLeft(Tau);

      FidLoss = this->CalculateFidelityLossLeft();
      LambdaDiff = this->CalculateLambdaDiffLeft();
   }

   if (Verbose > 0)
      std::cout << "FidelityLossLeft=" << FidLoss
                << " LambdaDiffLeft=" << LambdaDiff << std::endl;
}

void
IBC_TDVP::Evolve()
{
   ++TStep;
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();
   std::vector<double>::const_iterator GammaEnd = Comp.Gamma.cend();
   --GammaEnd;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowLeft();

   this->SweepLeftEW((*Gamma)*Timestep);
   ++Gamma;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowRight();

   while(Gamma != GammaEnd)
   {
      this->SweepRightEW((*Gamma)*Timestep);
      ++Gamma;

      this->SweepLeftEW((*Gamma)*Timestep);
      ++Gamma;
   }

   this->SweepRightFinalEW((*Gamma)*Timestep);
}

void
IBC_TDVP::EvolveExpand()
{
   ++TStep;
   Eps1SqSum = 0.0;
   Eps2SqSum = 0.0;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();
   std::vector<double>::const_iterator GammaEnd = Comp.Gamma.cend();
   --GammaEnd;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowLeft();

   this->SweepLeftExpandEW((*Gamma)*Timestep);
   ++Gamma;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowRight();

   while(Gamma != GammaEnd)
   {
      this->SweepRightEW((*Gamma)*Timestep);
      ++Gamma;

      this->SweepLeftExpandEW((*Gamma)*Timestep);
      ++Gamma;
   }

   this->SweepRightFinalEW((*Gamma)*Timestep);
}

void
IBC_TDVP::Evolve2()
{
   ++TStep;
   TruncErrSum = 0.0;

   std::vector<double>::const_iterator Gamma = Comp.Gamma.cbegin();

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowLeft();

   this->SweepLeft2((*Gamma)*Timestep);
   ++Gamma;

   if (NExpand != 0)
      if (TStep % NExpand == 0)
         this->ExpandWindowRight();

   this->SweepRight2((*Gamma)*Timestep);
   ++Gamma;

   while(Gamma != Comp.Gamma.cend())
   {
      this->SweepLeft2((*Gamma)*Timestep);
      ++Gamma;

      this->SweepRight2((*Gamma)*Timestep);
      ++Gamma;
   }

   this->CalculateEps();
}
