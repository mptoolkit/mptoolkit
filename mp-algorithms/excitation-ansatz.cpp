// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/excitation-ansatz.cpp
//
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp-algorithms/excitation-ansatz.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/transfer.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "tensor/tensor_eigen.h"
#include "wavefunction/operator_actions.h"

// The tolerance of the trace of the left/right boundary eigenvectors for
// fixing their relative phase.
double const TraceTol = 1e-8;

HEff::HEff(InfiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& HamMPO_,
           QuantumNumbers::QuantumNumber const& Q_, ProductMPO const& StringOp_,
           double k, double GMRESTol_, int Verbose_)
   : PsiLeft(Psi_), PsiRight(Psi_), HamMPO(HamMPO_), Q(Q_),
     StringOp(StringOp_), GMRESTol(GMRESTol_), Verbose(Verbose_)
{
   this->SetK(k);

   // Get PsiLeft and PsiRight as LinearWavefunctions.
   std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);

   MatrixOperator U;
   RealDiagonalOperator D;
   std::tie(U, D, PsiLinearRight) = get_right_canonical(PsiRight);
   PsiLinearRight.set_front(prod(U, PsiLinearRight.get_front()));

   // Get the leading eigenvectors for the mixed transfer matrix of PsiLeft
   // and PsiRight: for use with SolveSimpleMPOLeft/Right2.
   RhoL = delta_shift(PsiLeft.lambda_r(), PsiLeft.qshift());
   RhoR = PsiLeft.lambda_r();

   if (!StringOp.is_null())
   {
      // Calculate the left/right eigenvectors of the mixed transfer
      // matrices with the string operator corresponding to Ty.
      std::tie(std::ignore, TyLRLeft, TyLRRight) = get_transfer_eigenpair(PsiLinearLeft, PsiLinearRight, PsiLeft.qshift(), StringOp);
      TyLRRight = delta_shift(TyLRRight, PsiLeft.qshift());
      std::tie(std::ignore, TyRLLeft, TyRLRight) = get_transfer_eigenpair(PsiLinearRight, PsiLinearLeft, PsiLeft.qshift(), StringOp);
      TyRLLeft = delta_shift(TyRLLeft, adjoint(PsiRight.qshift()));
   }

   // Ensure HamMPO is the correct size.
   if (HamMPO.size() < PsiLeft.size())
      HamMPO = repeat(HamMPO, PsiLeft.size() / HamMPO.size());
   CHECK_EQUAL(HamMPO.size(), PsiLeft.size());

   // Solve the left Hamiltonian environment.
   BlockHamL = Initial_E(HamMPO, PsiLeft.Basis1());
   std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPO, GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Left energy = " << LeftEnergy << std::endl;

   BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

   // Solve the right Hamiltonian environment.
   BlockHamR = Initial_F(HamMPO, PsiLinearRight.Basis2());
   MatrixOperator Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
   Rho = delta_shift(Rho, adjoint(PsiRight.qshift()));

   std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiLinearRight, PsiRight.qshift(),
                                                           HamMPO, Rho, GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Right energy = " << RightEnergy << std::endl;

   // Remove the contribution from the ground state energy density.
   // To do this we need to remove the energy density contribution from one
   // unit cell due to the excitation ansatz window (the RightEnergy term),
   // and one contribution from the "bond energy", which is the energy
   // contribution from the terms in the Hamiltonian which cross the bond at
   // a unit cell boundary.
   std::complex<double> BondEnergy = inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamR);

   // An alternate way to calculate the bond energy using only the right
   // block Hamiltonian by essentially setting the upper-right element in
   // the unit cell MPO to be zero.
   // Unfortunately, this method does not work when there are diagonal terms
   // in the MPO.
#if 0
   std::vector<std::vector<int>> Mask = mask_row(HamMPO, 0);
   Mask.back().back() = 0;
   MatrixOperator C = inject_right_mask(BlockHamR, PsiLinearRight, HamMPO.data(), PsiLinearRight, Mask)[0];
   C.delta_shift(adjoint(PsiRight.qshift()));
   std::complex<double> BondEnergy = inner_prod(Rho, C);
#endif

   if (Verbose > 0)
      std::cout << "Bond energy = " << BondEnergy << std::endl;

   BlockHamR.front() -= (RightEnergy + BondEnergy) * BlockHamR.back();

   this->Initialize();
}

HEff::HEff(InfiniteWavefunctionLeft const& PsiLeft_, InfiniteWavefunctionLeft const& PsiRight_,
           BasicTriangularMPO const& HamMPO_, QuantumNumbers::QuantumNumber const& Q_,
           ProductMPO const& StringOp_, double k, double GMRESTol_, int Verbose_)
   : PsiLeft(PsiLeft_), PsiRight(PsiRight_), HamMPO(HamMPO_), Q(Q_),
     StringOp(StringOp_), GMRESTol(GMRESTol_), Verbose(Verbose_)
{
   CHECK_EQUAL(PsiLeft.size(), PsiRight.size());
   CHECK_EQUAL(PsiLeft.qshift(), PsiRight.qshift());

   this->SetK(k);

   // Get PsiLeft and PsiRight as LinearWavefunctions.
   std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);

   MatrixOperator U;
   RealDiagonalOperator D;
   std::tie(U, D, PsiLinearRight) = get_right_canonical(PsiRight);
   PsiLinearRight.set_front(prod(U, PsiLinearRight.get_front()));

   // Since the leading eigenvalue of the left/right mixed transfer matrix
   // has magnitude < 1, we do not need to orthogonalize the E/F matrix
   // elements against its eigenvectors.
   RhoL = MatrixOperator();
   RhoR = MatrixOperator();

   if (!StringOp.is_null())
   {
      TyLRLeft = MatrixOperator();
      TyLRRight = MatrixOperator();
      TyRLLeft = MatrixOperator();
      TyRLRight = MatrixOperator();
   }

   // Ensure HamMPO is the correct size.
   if (HamMPO.size() < PsiLeft.size())
      HamMPO = repeat(HamMPO, PsiLeft.size() / HamMPO.size());
   CHECK_EQUAL(HamMPO.size(), PsiLeft.size());

   // Solve the left Hamiltonian environment.
   BlockHamL = Initial_E(HamMPO, PsiLeft.Basis1());
   std::complex<double> LeftEnergy = SolveSimpleMPO_Left(BlockHamL, PsiLeft, HamMPO, GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Left energy = " << LeftEnergy << std::endl;

   BlockHamL = delta_shift(BlockHamL, adjoint(PsiLeft.qshift()));

   // Solve the right Hamiltonian environment.
   BlockHamR = Initial_F(HamMPO, PsiLinearRight.Basis2());
   MatrixOperator Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
   Rho = delta_shift(Rho, adjoint(PsiRight.qshift()));

   std::complex<double> RightEnergy = SolveSimpleMPO_Right(BlockHamR, PsiLinearRight, PsiRight.qshift(),
                                                           HamMPO, Rho, GMRESTol, Verbose-1);
   if (Verbose > 0)
      std::cout << "Right energy = " << RightEnergy << std::endl;

   // Solve the right Hamiltonian environment for PsiLeft to find the "bond
   // energy": see the comments above.
   LinearWavefunction PsiLinear;
   std::tie(U, D, PsiLinear) = get_right_canonical(PsiLeft);
   PsiLinear.set_front(prod(U, PsiLinear.get_front()));

   StateComponent BlockHamLR = Initial_F(HamMPO, PsiLinear.Basis2());
   Rho = scalar_prod(U*D*herm(U), herm(U*D*herm(U)));
   Rho = delta_shift(Rho, adjoint(PsiLeft.qshift()));

   SolveSimpleMPO_Right(BlockHamLR, PsiLinear, PsiLeft.qshift(), HamMPO, Rho, GMRESTol, Verbose-1);
   std::complex<double> BondEnergy = inner_prod(prod(PsiLeft.lambda_r(), prod(BlockHamL, PsiLeft.lambda_r())), BlockHamLR);

   if (Verbose > 0)
      std::cout << "Bond energy = " << BondEnergy << std::endl;

   // Remove the contribution from the ground state energy density.
   BlockHamR.front() -= (RightEnergy + BondEnergy) * BlockHamR.back();

   this->Initialize();
}

void
HEff::Initialize()
{
   // Get the null space matrices corresponding to each A-matrix in PsiLeft.
   for (StateComponent C : PsiLinearLeft)
      NullLeftDeque.push_back(NullSpace2(C));

   // Construct the partially contracted left Hamiltonian environments in the unit cell.
   BlockHamLDeque.push_back(BlockHamL);
   auto CL = PsiLinearLeft.begin();
   auto O = HamMPO.begin();
   while (CL != PsiLinearLeft.end())
   {
      BlockHamLDeque.push_back(contract_from_left(*O, herm(*CL), BlockHamLDeque.back(), *CL));
      ++CL, ++O;
   }

   // Same for the right environments.
   BlockHamRDeque.push_front(BlockHamR);
   auto CR = PsiLinearRight.end();
   O = HamMPO.end();
   while (CR != PsiLinearRight.begin())
   {
      --CR, --O;
      BlockHamRDeque.push_front(contract_from_right(herm(*O), *CR, BlockHamRDeque.front(), herm(*CR)));
   }

   if (!StringOp.is_null())
   {
      MatrixOperator Tmp;

      // Calculate the left/right eigenvectors of the transfer matrix with
      // the string operator corresponding to Ty, fixing their norms and
      // relative phases.
      std::tie(std::ignore, TyL, Tmp) = get_transfer_eigenpair(PsiLinearLeft, PsiLinearLeft, PsiLeft.qshift(), StringOp);

      // Normalize TyL s.t. the sum of the singular values of Tmp = 1.
      MatrixOperator U, Vh;
      RealDiagonalOperator D;
      SingularValueDecomposition(Tmp, U, D, Vh);

      Tmp *= 1.0 / trace(D);
      TyL *= 1.0 / inner_prod(delta_shift(Tmp, PsiLeft.qshift()), TyL);

      std::tie(std::ignore, Tmp, TyR) = get_transfer_eigenpair(PsiLinearRight, PsiLinearRight, PsiRight.qshift(), StringOp);

      // Normalize TyR s.t. the sum of the singular values of Tmp = 1.
      SingularValueDecomposition(Tmp, U, D, Vh);

      Tmp *= 1.0 / trace(D);
      TyR *= 1.0 / inner_prod(delta_shift(TyR, PsiRight.qshift()), Tmp);

      // Fix the phases of TyL and TyR by setting the phases of their traces to be zero.
      // TODO: Figure out a better method to fix the phase: see wavefunction/ibc.cpp.
      if (TyL.Basis1() == TyL.Basis2() && TyR.Basis1() == TyR.Basis2())
      {
         std::complex<double> TyLTrace = trace(TyL);

         if (std::abs(TyLTrace) > TraceTol)
            TyL *= std::conj(TyLTrace) / std::abs(TyLTrace);
         else
            WARNING("The trace of TyL is below threshold, so Ty will have a spurious phase contribution.")(TyLTrace);

         std::complex<double> TyRTrace = trace(TyR);

         if (std::abs(TyRTrace) > TraceTol)
            TyR *= std::conj(TyRTrace) / std::abs(TyRTrace);
         else
            WARNING("The trace of TyR is below threshold, so Ty will have a spurious phase contribution.")(TyRTrace);
      }
      else
         WARNING("Psi1 and Psi2 have different boundary bases, so Ty will have a spurious phase contribution.");

      // Only needed when adding TyEff to HEff.
#if 0
      // Construct the partially contracted versions of TyL and TyR.
      StateComponent TyLSC = StateComponent(StringOp.Basis1(), PsiLeft.Basis1(), PsiLeft.Basis1());
      TyLSC.front() = TyL;

      TyLDeque.push_back(TyLSC);

      CL = PsiLinearLeft.begin();
      O = StringOp.begin();
      while (CL != PsiLinearLeft.end())
      {
         TyLDeque.push_back(contract_from_left(*O, herm(*CL), TyLDeque.back(), *CL));
         ++CL, ++O;
      }

      StateComponent TyRSC = StateComponent(StringOp.Basis2(), PsiRight.Basis2(), PsiRight.Basis2());
      TyRSC.front() = TyR;

      TyRDeque.push_front(TyRSC);
      CR = PsiLinearRight.end();
      O = StringOp.end();
      while (CR != PsiLinearRight.begin())
      {
         --CR, --O;
         TyRDeque.push_front(contract_from_right(herm(*O), *CR, TyRDeque.front(), herm(*CR)));
      }
#endif
   }
}

std::deque<MatrixOperator>
HEff::operator()(std::deque<MatrixOperator> const& XDeque) const
{
   std::deque<StateComponent> BDeque = this->ConstructBDeque(XDeque);
   LinearWavefunction PsiTri = this->ConstructPsiTri(BDeque);

   // Calcaulate the terms in the triangular E and F matrices where there is
   // one B-matrix on the top.
   StateComponent BlockHamLTri, BlockHamRTri;

   SolveSimpleMPO_Left2(BlockHamLTri, BlockHamL, PsiLinearLeft, PsiLinearRight, PsiTri,
                        PsiLeft.qshift(), HamMPO, RhoL, RhoL, ExpIK, GMRESTol, Verbose-1);
   SolveSimpleMPO_Right2(BlockHamRTri, BlockHamR, PsiLinearLeft, PsiLinearRight, PsiTri,
                         PsiRight.qshift(), HamMPO, RhoR, RhoR, ExpIK, GMRESTol, Verbose-1);

   // Shift the phases by one unit cell.
   BlockHamLTri *= ExpIK;
   BlockHamRTri *= ExpIK;

   // Calculate the contribution to HEff corresponding to where the
   // B-matrices are on the same site.
   std::deque<MatrixOperator> Result;
   auto B = BDeque.begin();
   auto NL = NullLeftDeque.begin();
   auto O = HamMPO.begin();
   auto BHL = BlockHamLDeque.begin();
   auto BHR = BlockHamRDeque.begin();
   ++BHR;
   while (B != BDeque.end())
   {
      Result.push_back(scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *NL)), *BHR));
      ++B, ++NL, ++O, ++BHL, ++BHR;
   }

   // Calculate the contribution where the top B-matrix is in the left
   // semi-infinite part.
   StateComponent Tmp = BlockHamLTri;
   B = BDeque.begin();
   NL = NullLeftDeque.begin();
   auto CL = PsiLinearLeft.begin();
   auto CR = PsiLinearRight.begin();
   O = HamMPO.begin();
   BHL = BlockHamLDeque.begin();
   BHR = BlockHamRDeque.begin();
   ++BHR;
   auto R = Result.begin();
   while (B != BDeque.end())
   {
      *R += scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *NL)), *BHR);
      ++R;
      if (R != Result.end())
         Tmp = contract_from_left(*O, herm(*CR), Tmp, *CL) + contract_from_left(*O, herm(*B), *BHL, *CL);
      ++B, ++NL, ++CL, ++CR, ++O, ++BHL, ++BHR;
   }

   // Calculate the contribution where the top B-matrix is in the right
   // semi-infinite part.
   Tmp = BlockHamRTri;
   B = BDeque.end();
   NL = NullLeftDeque.end();
   CL = PsiLinearLeft.end();
   CR = PsiLinearRight.end();
   O = HamMPO.end();
   BHL = BlockHamLDeque.end();
   --BHL;
   BHR = BlockHamRDeque.end();
   R = Result.end();
   while (B != BDeque.begin())
   {
      --B, --NL, --CL, --CR, --O, --BHL, --BHR, --R;
      *R += scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *NL)), Tmp);
      if (B != BDeque.begin())
         Tmp = contract_from_right(herm(*O), *CL, Tmp, herm(*CR)) + contract_from_right(herm(*O), *B, *BHR, herm(*CR));
   }

   // Code to target excitations with a specific y-momentum (WIP).
#if 0
   if (!StringOp.is_null())
   {
      std::complex<double> Alpha = 10.0;

      MatrixOperator E, F;
      SolveStringMPO_Left2(E, TyL, PsiLinearLeft, PsiLinearRight, PsiTri,
                           PsiLeft.qshift(), StringOp, TyLRLeft, TyLRRight, ExpIK, GMRESTol, Verbose-1);
      SolveStringMPO_Right2(F, TyR, PsiLinearLeft, PsiLinearRight, PsiTri,
                            PsiRight.qshift(), StringOp, TyRLLeft, TyRLRight, ExpIK, GMRESTol, Verbose-1);

      E *= ExpIK;
      F *= ExpIK;

      // Calculate the contribution to HEff corresponding to where the
      // B-matrices are on the same site.
      B = BDeque.begin();
      NL = NullLeftDeque.begin();
      O = StringOp.begin();
      BHL = TyLDeque.begin();
      BHR = TyRDeque.begin();
      ++BHR;
      R = Result.begin();
      while (B != BDeque.end())
      {
         *R += -Alpha * scalar_prod(herm(contract_from_left(*O, herm(*B), *BHL, *NL)), *BHR);
         ++B, ++NL, ++O, ++BHL, ++BHR, ++R;
      }

      // Calculate the contribution where the top B-matrix is in the left
      // semi-infinite part.
      Tmp = StateComponent(StringOp.Basis1(), PsiLinearRight.Basis1(), PsiLinearLeft.Basis1());
      Tmp.front() = E;
      B = BDeque.begin();
      NL = NullLeftDeque.begin();
      CL = PsiLinearLeft.begin();
      CR = PsiLinearRight.begin();
      O = StringOp.begin();
      BHL = TyLDeque.begin();
      BHR = TyRDeque.begin();
      ++BHR;
      R = Result.begin();
      while (B != BDeque.end())
      {
         *R += -Alpha * scalar_prod(herm(contract_from_left(*O, herm(*CR), Tmp, *NL)), *BHR);
         ++R;
         if (R != Result.end())
            Tmp = contract_from_left(*O, herm(*CR), Tmp, *CL) + contract_from_left(*O, herm(*B), *BHL, *CL);
         ++B, ++NL, ++CL, ++CR, ++O, ++BHL, ++BHR;
      }

      // Calculate the contribution where the top B-matrix is in the right
      // semi-infinite part.
      Tmp = StateComponent(StringOp.Basis2(), PsiLinearLeft.Basis2(), PsiLinearRight.Basis2());
      Tmp.front() = F;
      B = BDeque.end();
      NL = NullLeftDeque.end();
      CL = PsiLinearLeft.end();
      CR = PsiLinearRight.end();
      O = StringOp.end();
      BHL = TyLDeque.end();
      --BHL;
      BHR = TyRDeque.end();
      R = Result.end();
      while (B != BDeque.begin())
      {
         --B, --NL, --CL, --CR, --O, --BHL, --BHR, --R;
         *R += -Alpha * scalar_prod(herm(contract_from_left(*O, herm(*CL), *BHL, *NL)), Tmp);
         if (B != BDeque.end())
            Tmp = contract_from_right(herm(*O), *CL, Tmp, herm(*CR)) + contract_from_right(herm(*O), *B, *BHR, herm(*CR));
      }
   }
#endif

   return Result;
}

std::complex<double>
HEff::Ty(std::deque<MatrixOperator> const& XDeque) const
{
   std::deque<StateComponent> BDeque = this->ConstructBDeque(XDeque);
   LinearWavefunction PsiTri = this->ConstructPsiTri(BDeque);

   MatrixOperator E, F;
   SolveStringMPO_Left2(E, TyL, PsiLinearLeft, PsiLinearRight, PsiTri,
                        PsiLeft.qshift(), StringOp, TyLRLeft, TyLRRight, ExpIK, GMRESTol, Verbose-1);
   SolveStringMPO_Right2(F, TyR, PsiLinearLeft, PsiLinearRight, PsiTri,
                         PsiRight.qshift(), StringOp, TyRLLeft, TyRLRight, ExpIK, GMRESTol, Verbose-1);

   E *= ExpIK;
   F *= ExpIK;

   std::complex<double> Ty = inner_prod(inject_left(TyL, PsiTri, StringOp, PsiTri), TyR)
                           + inner_prod(inject_left(E, PsiLinearRight, StringOp, PsiTri), TyR)
                           + inner_prod(inject_left(TyL, PsiLinearLeft, StringOp, PsiTri), F);
   return Ty;
}

std::deque<StateComponent>
HEff::ConstructBDeque(std::deque<MatrixOperator> const& XDeque) const
{
   std::deque<StateComponent> BDeque;
   auto NL = NullLeftDeque.begin();
   auto X = XDeque.begin();
   while (NL != NullLeftDeque.end())
   {
      BDeque.push_back(prod(*NL, *X));
      ++NL, ++X;
   }

   return BDeque;
}

LinearWavefunction
HEff::ConstructPsiTri(std::deque<StateComponent> const& BDeque) const
{
   LinearWavefunction PsiTri;

   if (PsiLeft.size() == 1)
      PsiTri.push_back(BDeque.back());
   else
   {
      auto CL = PsiLinearLeft.begin();
      auto CR = PsiLinearRight.begin();
      auto B = BDeque.begin();
      SumBasis<VectorBasis> NewBasis0((*CL).Basis2(), (*B).Basis2());
      PsiTri.push_back(tensor_row_sum(*CL, *B, NewBasis0));
      ++CL, ++CR, ++B;
      for (int i = 1; i < PsiLeft.size()-1; ++i)
      {
         StateComponent Z = StateComponent((*CL).LocalBasis(), (*CR).Basis1(), (*CL).Basis2());
         SumBasis<VectorBasis> NewBasis1((*CL).Basis2(), (*B).Basis2());
         SumBasis<VectorBasis> NewBasis2((*CL).Basis1(), (*CR).Basis1());
         PsiTri.push_back(tensor_col_sum(tensor_row_sum(*CL, *B, NewBasis1), tensor_row_sum(Z, *CR, NewBasis1), NewBasis2));
         ++CL, ++CR, ++B;
      }
      SumBasis<VectorBasis> NewBasis3((*B).Basis1(), (*CR).Basis1());
      PsiTri.push_back(tensor_col_sum(*B, *CR, NewBasis3));
   }

   return PsiTri;
}

std::deque<MatrixOperator>
HEff::InitialGuess() const
{
   std::deque<MatrixOperator> Result;
   auto NL = NullLeftDeque.begin();
   auto CR = PsiLinearRight.begin();
   while (NL != NullLeftDeque.end())
   {
      MatrixOperator C = MakeRandomMatrixOperator((*NL).Basis2(), (*CR).Basis2(), Q);
      C *= 1.0 / norm_frob(C);
      Result.push_back(C);
      ++NL, ++CR;
   }

   return Result;
}

std::deque<PackMatrixOperator>
HEff::PackInitialize() const
{
   std::deque<PackMatrixOperator> Result;
   auto NL = NullLeftDeque.begin();
   auto CR = PsiLinearRight.begin();
   while (NL != NullLeftDeque.end())
   {
      Result.push_back(PackMatrixOperator((*NL).Basis2(), (*CR).Basis2(), Q));
      ++NL, ++CR;
   }

   return Result;
}

EAWavefunction
HEff::ConstructEAWavefunction(std::deque<MatrixOperator> XDeque) const
{
   std::vector<WavefunctionSectionLeft> WindowVec;
   auto NL = NullLeftDeque.begin();
   auto X = XDeque.begin();
   while (NL != NullLeftDeque.end())
   {
      LinearWavefunction Psi;
      Psi.push_back(*NL);
      WindowVec.push_back(WavefunctionSectionLeft::ConstructFromLeftOrthogonal(std::move(Psi), *X, Verbose-1));
      ++NL, ++X;
   }

   // TODO: Do we need to update PsiLeft/PsiRight?
   return EAWavefunction(PsiLeft, WindowVec, PsiRight, ExpIK);
}

void
HEff::SetK(double k)
{
   ExpIK = exp(std::complex<double>(0.0, math_const::pi) * k);
}


PackHEff::PackHEff(HEff H_)
    : H(H_)
{
   Pack = H.PackInitialize();
   Size = 0;
   for (auto P : Pack)
      Size += P.size();
}

void
PackHEff::operator()(std::complex<double> const* In_, std::complex<double>* Out_) const
{
   std::deque<MatrixOperator> XDeque = this->unpack(In_);
   XDeque = H(XDeque);
   this->pack(XDeque, Out_);
}

std::deque<MatrixOperator>
PackHEff::unpack(std::complex<double> const* In_) const
{
   std::complex<double> const* In = In_;
   std::deque<MatrixOperator> XDeque;
   for (auto P : Pack)
   {
      XDeque.push_back(P.unpack(In));
      In += P.size();
   }

   return XDeque;
}

void
PackHEff::pack(std::deque<MatrixOperator> XDeque, std::complex<double>* Out_) const
{
   std::complex<double>* Out = Out_;
   auto P = Pack.begin();
   auto X = XDeque.begin();
   while (P != Pack.end())
   {
      (*P).pack(*X, Out);
      Out += (*P).size();
      ++P, ++X;
   }
}
