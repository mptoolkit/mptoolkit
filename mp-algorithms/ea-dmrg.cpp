// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ea-dmrg.cpp
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "ea-dmrg.h"
#include "linearalgebra/arpack_wrapper.h"

// Tolerance for checking whether a window is orthogonal to the ground state.
double const OrthoTol = 1e-4;

EA_DMRG::EA_DMRG(EAWavefunction const& Psi_, BasicTriangularMPO const& HamMPO, EA_DMRGSettings Settings_)
   : Psi(Psi_), Tol(Settings_.Tol), Quiet(Settings_.Quiet), Verbose(Settings_.Verbose)
{
   PRECONDITION(Psi.window_size() > 1);

   // Extract the left and right semi-infinite boundaries.
   InfiniteWavefunctionLeft PsiLeft = Psi.left();
   inplace_qshift(PsiLeft, Psi.left_qshift());
   PsiLeft.rotate_left(Psi.left_index());

   InfiniteWavefunctionRight PsiRight = Psi.right();
   inplace_qshift(PsiRight, Psi.right_qshift());
   PsiRight.rotate_left(Psi.right_index());

   // Extract the windows.
   for (WavefunctionSectionLeft Window : Psi.window_vec())
   {
      LinearWavefunction PsiLinear;
      MatrixOperator U;
      std::tie(PsiLinear, U) = get_left_canonical(Window);
      PsiLinear.set_back(PsiLinear.get_back()*U);
      WindowVec.push_back(PsiLinear);
   }

   // Make sure the windows obey the left-gauge fixing condition.
   LinearWavefunction PsiLinearLeft;
   std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);

   auto Window = WindowVec.begin();
   auto C = PsiLinearLeft.begin();
   while (Window != WindowVec.end())
   {
      StateComponent WFront = Window->get_front();
      Window->pop_front();

      StateComponent N = NullSpace2(*C);

      // Make sure the first site each of window is of the form NX.
      if (norm_frob(scalar_prod(herm(*C), WFront)) > OrthoTol)
         throw std::runtime_error("EA_DMRG: fatal: window does not obey the left-gauge fixing condition.");

      // Absorb X into the second site.
      MatrixOperator X = scalar_prod(herm(N), WFront);
      Window->set_front(prod(X, Window->get_front()));

      Window->push_front(N);

      WIVec.push_back(Window->end());
      --WIVec.back();

      ++Window, ++C;
   }

   int WS = Psi.window_size();

   LeftStop = 1;
   RightStop = WS-1;
   Site = RightStop;

   // Initialize the EFMatrix class.
   EFMatrixSettings Settings;
   Settings.Tol = Settings_.GMRESTol;
   Settings.UnityEpsilon = Settings_.UnityEpsilon;
   Settings.EAOptimization = true;
   Settings.SubtractEnergy = true;
   Settings.Verbose = Verbose;

   EF = EFMatrix(HamMPO, Settings);
   EF.SetPsi(false, PsiLeft);
   EF.SetPsi(true, PsiRight, Psi.exp_ik());
   EF.SetWindowUpper(WindowVec);
   EF.SetWindowLower(WindowVec);

   // FIXME: If we do not run GetHEff before using the ARPACK wrapper, we run into memory issues.
   std::deque<StateComponent> HWDeque = EF.GetHEff(Site);
   // Calculate and print the initial energy.
   if (!Quiet)
   {
      double InitialE = 0.0;
      auto WI = WIVec.begin();
      auto HW = HWDeque.begin();
      while (HW != HWDeque.end())
      {
         InitialE += std::real(inner_prod(*(*WI), *HW));
         ++WI, ++HW;
      }
      std::cout << "Sweep=" << Sweep
                << " InitialE=" << InitialE << std::endl;
   }
}

EAWavefunction
EA_DMRG::Wavefunction() const
{
   // Convert the windows to WavefunctionSectionLefts.
   std::vector<WavefunctionSectionLeft> WSLVec;
   for (auto Window : WindowVec)
   {
      MatrixOperator Identity = MatrixOperator::make_identity(Window.Basis2());
      WSLVec.push_back(WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Window), Identity, Verbose));
   }

   EAWavefunction Result(Psi.left(), WSLVec, Psi.right(), Psi.left_qshift(), Psi.right_qshift(),
                         Psi.left_index(), Psi.right_index(), Psi.exp_ik(), Psi.gs_overlap());

   // Stream the boundaries if the initial state does.
   Result.set_left_filename(Psi.get_left_filename());
   Result.set_right_filename(Psi.get_right_filename());

   return Result;
}

void
EA_DMRG::SolveCurrentSite()
{
   // Get the window components for this window site.
   std::deque<StateComponent> WDeque;
   for (auto const& WI : WIVec)
      WDeque.push_back(*WI);

   // Set the effective Hamiltonian functors.
   HEff H(&EF, Site);
   PackHEff PackH(&H, WDeque);

   // Pack the initial window components.
   std::vector<std::complex<double>> Init(PackH.size());
   PackH.pack(WDeque, Init.data());

   // Solve for the minimum of HEff using ARPACK.
   std::vector<std::complex<double>> Solution;
   LinearAlgebra::Vector<std::complex<double>> EValues
      = LinearAlgebra::DiagonalizeARPACK(PackH, PackH.size(), 1, LinearAlgebra::WhichEigenvalues::SmallestReal,
                                         Init.data(), Tol, &Solution, 0, true, Verbose);

   if (!Quiet)
      std::cout << "Sweep=" << Sweep
                << " Site=" << Site
                << " E=" << std::real(EValues[0])
                << std::endl;

   // Unpack the solution.
   WDeque = PackH.unpack(Solution.data());

   // Update the windows.
   auto WI = WIVec.begin();
   auto W = WDeque.begin();
   while (W != WDeque.end())
   {
      *(*WI) = *W;
      ++WI, ++W;
   }
}

void
EA_DMRG::IterateLeft()
{
   CHECK(Site > LeftStop);

   std::deque<StateComponent> WDeque;
   for (auto& WI : WIVec)
   {
      // Perform SVD to right-orthogonalize current site.
      MatrixOperator U;
      RealDiagonalOperator D;

      std::tie(U, D) = OrthogonalizeBasis1(*WI);
      WDeque.push_back(*WI);

      // Move to the next site.
      --WI;

      *WI = prod(*WI, U*D);
   }

   // Update the windows.
   EF.SetWindowUpper(Site+1, WDeque);
   EF.SetWindowLower(Site+1, WDeque);

   --Site;
}

void
EA_DMRG::IterateRight()
{
   CHECK(Site < RightStop);

   std::deque<StateComponent> WDeque;
   for (auto& WI : WIVec)
   {
      // Perform SVD to left-orthogonalize current site.
      MatrixOperator Vh;
      RealDiagonalOperator D;

      std::tie(D, Vh) = OrthogonalizeBasis2(*WI);
      WDeque.push_back(*WI);

      // Move to the next site.
      ++WI;

      *WI = prod(D*Vh, *WI);
   }

   // Update the windows.
   EF.SetWindowUpper(Site+1, WDeque);
   EF.SetWindowLower(Site+1, WDeque);

   ++Site;
}

void
EA_DMRG::SweepLR()
{
   CHECK(Site == RightStop);
   ++Sweep;

   while (Site > LeftStop)
   {
      this->IterateLeft();
      this->SolveCurrentSite();
   }

   ++Sweep;

   while (Site < RightStop)
   {
      this->IterateRight();
      this->SolveCurrentSite();
   }
}

HEff::HEff(EFMatrix* EF_, int Site_)
   : EF(EF_), Site(Site_)
{
}

std::deque<StateComponent>
HEff::operator()(std::deque<StateComponent> WDeque)
{
   EF->SetWindowLower(Site+1, WDeque);
   return EF->GetHEff(Site);
}

PackHEff::PackHEff(HEff* H_, std::deque<StateComponent> WDeque)
   : H(H_)
{
   for (auto const& W : WDeque)
      Pack.push_back(PackStateComponent(W));
   Size = 0;
   for (auto P : Pack)
      Size += P.size();
}

void
PackHEff::operator()(std::complex<double> const* In_, std::complex<double>* Out_) const
{
   std::deque<StateComponent> WDeque = this->unpack(In_);
   WDeque = (*H)(WDeque);
   this->pack(WDeque, Out_);
}

std::deque<StateComponent>
PackHEff::unpack(std::complex<double> const* In_) const
{
   std::complex<double> const* In = In_;
   std::deque<StateComponent> WDeque;
   for (auto P : Pack)
   {
      WDeque.push_back(P.unpack(In));
      In += P.size();
   }
   return WDeque;
}

void
PackHEff::pack(std::deque<StateComponent> WDeque, std::complex<double>* Out_) const
{
   std::complex<double>* Out = Out_;
   auto P = Pack.begin();
   auto W = WDeque.begin();
   while (P != Pack.end())
   {
      std::complex<double>* Tmp = (*P).pack(*W, Out);
      Out += (*P).size();
      ++P, ++W;
   }
}
