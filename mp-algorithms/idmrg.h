// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/finitedmrg.h
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
//
// iDMRG with environment expansion
//
// The wavefunction has a unit cell of N sites. To allow for pre-expansion, the window
// needs to include one additional site at each end of the lattice, so the window size is N+2.

// From an initial window of size N, can we bootstrap it?
// Start from
// E EEH F
//   AAC
// We can construct the final E-matrix, which gives the initial saved E, and saved A. If the unit cell is 1-site, then
// it will naturally end up as a copy (with delta-shift).
// EE EEH F
//  A AAC
// Sweep to the left. If we get to site 0 and we are at the beginning of the window, then we need to extend the
// window by one site. We get to
// EE HFF F
//  A CBB
// and we can do pre-expansion no problem. We end up with
// EE  HFF F
//  A λBBB
//
// Aside: at this point, we could solve the fixed-point equations
// AAA
// BBB
// from the previous sweep to the current sweep, and update the F-matrix environment.
//
// We can now obtain the final F matrix, to give
// EE  FFF F
//  A λBBB
// at which point we can update the saved right hand matrices, FF and λB
// I propose that we save F and λ together, as the usual saved environment, and the second F and B as the pre-expansion
//

//
// To initialize the system, we find the fixed point Hamiltonian for the left and right,
// and start from the state, inn a left-orthogonalized unit cell to perform a right-to-left sweep:
//
// EE EEH F
//  A AAC B
// where EE is the two sites beyond the left edge of the unit cell, F is the site beyond the right edge,
// AAC is the wavefunction, A is left-orthogonalized, C is the orthogonality center.  EEH is the
// environment for the unit cell, EE are E-matrices, H represents the Hamiltonian MPO at the
// orthogonality center.  We also need the Hamiltonian MPO for the unit cell plus one buffer site on each end.
// Note that the initial F environment only needs to be one site, since F is never re-used in the basis expansion.
// We also stash EE, A_E, λ_E as the saved left environment, and λ_F as the saved right environment.
//
// We then do conventional DMRG, with pre- and/or post-expansion, until we get to the left hand edge of the unit cell,
// giving
// EE HFF F
//    CBB
// At this point we update the left environment, replacing the EE matrices with the saved left environment,
// update C via C' = λ_E (λ_F)⁻¹ C
// Now we do the iteration as usual, which could involve pre-expansion (or even a 2-site update).
// Finally we do the truncation (possibly with post-expansion), to obtain
// EE  HFF F
//  A λBBB
// where we have done truncation and post-expansion to obtain λB.  We then calculate the final F matrix, and save
// λ as λ_F, B_F and FF as the saved right environment. We don't use F and λ_F immediately, in principle this could be done
// in another thread? Note that Basis1 of λ_F is the basis of A_E, and Basis2() of λ_F is B_F.
//
// We now sweep to the right, starting from
// EE HFF FF
//  A CBB B
// until we get to
// EE HFF F
//  A AAC B
// at which point we update the F matrices to λ_F B_F and FF via C' = C (λ_E)⁻¹ λ_F,
// update the right side using our saved λ_F and FF, as C' = C  which gives
// EE HFF FF
//  A AAC B
// We can now sweep to the left, and simultaneously save the left block, as C = A_E λ_E, and the new EE
//
// To allow pre-expansion, we need an additional environment matrix and corresponding A-matrix.

//
// In this version, we don't bother with pre-expansion.


#if !defined(MPTOOLKIT_MP_ALGORITHMS_IDMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_IDMRG_H

#include "mp-algorithms/dmrg.h"
#include "wavefunction/infinitewavefunctionleft.h"

class iDMRG : public DMRG
{
   public:
      iDMRG(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_);

      void UpdateLeftBlock();

      void UpdateRightBlock();

      // get the current wavefunction. This can be called at any time, and the wavefunction will be
      // converted to canonical form
      InfiniteWavefunctionLeft Wavefunction() const;

      virtual void check_structure() const;

      LinearWavefunction::const_iterator FirstSite, LastSite;

      QuantumNumber QShift;

      StateComponent SaveLeftHamiltonian;
      RealDiagonalOperator SaveLambda2;
      MatrixOperator SaveU2;

      StateComponent SaveRightHamiltonian;
      MatrixOperator SaveU1;
      RealDiagonalOperator SaveLambda1;
};


iDMRG::iDMRG(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_)
   : DMRG(Verbose_)
{
   LinearWavefunction Psi;
   RealDiagonalOperator Lambda;
   std::tie(Psi, Lambda) = get_left_canonical(Psi_);

   // Find the fixed point matrix elements
   StateComponent BlockHamL = Initial_E(Ham_, Psi.Basis1());
   std::cout << "Solving fixed-point Hamiltonian..." << std::endl;
   MatrixOperator Rho = delta_shift(scalar_prod(Lambda, herm(Lambda)), QShift);
   InitialEnergy = SolveHamiltonianMPO_Left(BlockHamL, Psi, QShift, HamMPO,
                                       Rho, GMRESTol, Verbose);
   std::cout << "Starting energy (left eigenvalue) = " << InitialEnergy << std::endl;


   StateComponent BlockHamR = Initial_F(HamMPO, Psi.Basis2());
   LinearWavefunction PsiR;
   RealDiagonalOperator D;
   std::tie(D, PsiR) = get_right_canonical(Psi_);

   BlockHamR = Initial_F(HamMPO, PsiR.Basis2());

   Rho = D*D;

   // We obtained Rho from the left side, so we need to delta shift to the right basis
   Rho = delta_shift(Rho, adjoint(QShift));

   std::complex<double> Energy = SolveHamiltonianMPO_Right(BlockHamR, PsiR, QShift, Ham_,
                                                           Rho, GMRESTol, Verbose);
   std::cout << "Starting energy (right eigenvalue) = " << Energy << std::endl;

   DMRG::Initialize(Psi, Ham_, BlockHamL, BlockHamR);

   // The InitialEnergy is the energy per unit cell.  At this point,
   // the Hamiltonian we have represents the sum of interactions across
   // one unit cell, but the boundary interactions are counted at both boundaries.
   // So this gives an energy that is (UnitCellSize+1) times the energy per site.
   // To compensate, we subtract a term off the RightHamiltonian
   // so that the initial energy is correct.


}

void
iDMRG::UpdateLeftBlock(double HMix)
{
   if (Verbose > 2)
   {
      std::cerr << "Updating left Hamiltonian.  Old basis has "
                << LeftHamiltonian.back().Basis1().total_dimension()
                << " states, new basis has states = "
                << SaveLeftHamiltonian.Basis1().total_dimension() << "\n";
   }

   StateComponent E = LeftHamiltonian.back();
   LeftHamiltonian = std::deque<StateComponent>(1, delta_shift(SaveLeftHamiltonian, QShift));

   MatrixOperator T = Solve_D_U_DInv(delta_shift(SaveLambda2, QShift), delta_shift(SaveU2, QShift), SaveLambda1);

   T = T * herm(SaveU1);

   (*C) = prod(T, *C);

   // normalize
   *C *= 1.0 / norm_frob(*C);

   // Subtract off the energy
   LeftHamiltonian.back().back() -= Solver_.LastEnergy() * LeftHamiltonian.back().front();

   this->CheckConsistency();
}

void
iDMRG::SaveLeftBlock(StatesInfo const& States)
{
   //TRACE(LeftHamiltonian.back());
   // When we save the block, we need to end up with
   // C.Basis2() == SaveLeftHamiltonian.Basis()
   CHECK(C == LastSite);
   StateComponent L = *C;
   std::tie(SaveLambda2, SaveU2) = SubspaceExpandBasis2(L, *H, LeftHamiltonian.back(),
                                                          MixingInfo, States, Info, RightHamiltonian.front());

   if (Verbose > 1)
   {
      std::cerr << "Saving left block for idmrg, states=" << Info.KeptStates()
                << " L.Basis2() is " << L.Basis2().total_dimension() << '\n';
   }
   SaveLeftHamiltonian = contract_from_left(*H, herm(L), LeftHamiltonian.back(), L);
   this->CheckConsistency();

   //TRACE(SaveLeftHamiltonian);
}



#endif
