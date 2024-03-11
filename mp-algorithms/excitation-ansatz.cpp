// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/excitation-ansatz.cpp
//
// Copyright (C) 2022-2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "mp-algorithms/excitation-ansatz.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/transfer.h"
#include "mp-algorithms/triangular_mpo_solver.h"
#include "tensor/tensor_eigen.h"
#include "wavefunction/operator_actions.h"

HEff::HEff(InfiniteWavefunctionLeft const& PsiLeft_, InfiniteWavefunctionRight const& PsiRight_,
           BasicTriangularMPO const& HamMPO_, EASettings const& Settings_)
   : PsiLeft(PsiLeft_), PsiRight(PsiRight_), HamMPO(HamMPO_),
     StringOp(Settings_.StringOp), GMRESTol(Settings_.GMRESTol),
     UnityEpsilon(Settings_.UnityEpsilon), kyFix(Settings_.kyFix),
     Alpha(Settings_.Alpha), Verbose(Settings_.Verbose)
{
   CHECK_EQUAL(PsiLeft.size(), PsiRight.size());
   CHECK_EQUAL(PsiLeft.qshift(), PsiRight.qshift());

   //this->SetK(Settings_.k); // We cannot use this method as it modifies EF, which hasn't been intialized yet.
   ExpIK = std::exp(std::complex<double>(0.0, math_const::pi) * Settings_.k);
   this->SetKY(Settings_.ky);

   // Ensure HamMPO is the correct size.
   if (HamMPO.size() < PsiLeft.size())
      HamMPO = repeat(HamMPO, PsiLeft.size() / HamMPO.size());
   CHECK_EQUAL(HamMPO.size(), PsiLeft.size());

   EFMatrixSettings Settings;
   Settings.Tol = GMRESTol;
   Settings.UnityEpsilon = UnityEpsilon;
   Settings.EAOptimization = true;
   Settings.SubtractEnergy = true;
   Settings.PhaseWarnings = Settings_.PhaseWarnings;
   Settings.Verbose = Verbose;

   EF = EFMatrix(HamMPO, Settings);
   EF.SetPsi({0}, PsiLeft);
   EF.SetPsi({Infinity}, PsiRight, ExpIK);
   EF.CalculateAllTEVs();

   if (!StringOp.is_null())
   {
      // Ensure StringOp is the correct size.
      if (StringOp.size() < PsiLeft.size())
         StringOp = repeat(StringOp, PsiLeft.size() / StringOp.size());
      CHECK_EQUAL(StringOp.size(), PsiLeft.size());

      EFMatrixSettings SettingsTy;
      SettingsTy.Tol = GMRESTol;
      SettingsTy.UnityEpsilon = UnityEpsilon;
      SettingsTy.PhaseWarnings = Settings_.PhaseWarnings;
      SettingsTy.Verbose = Verbose;

      EFTy = EFMatrix(StringOp, SettingsTy);
      EFTy.SetPsi({0}, PsiLeft);
      EFTy.SetPsi({Infinity}, PsiRight, ExpIK);
      EFTy.CalculateAllTEVs();
   }

   // Get PsiLeft and PsiRight as LinearWavefunctions.
   std::tie(PsiLinearLeft, std::ignore) = get_left_canonical(PsiLeft);
   std::tie(std::ignore, PsiLinearRight) = get_right_canonical(PsiRight);

   // Check that there are compatible quantum number sectors between PsiLeft
   // and PsiRight at the unit cell boundary.
   {
      auto Test = PackMatrixOperator(MatrixOperator(PsiLinearLeft.Basis1(), PsiLinearRight.Basis1()));
      if (Test.size() == 0)
      {
         std::string ErrorMessage = "fatal: The effective Hamiltonian has dimension zero. "
                                    "This probably means the bases of the left and right wavefunctions have incompatible quantum number sectors: "
                                    // This error message may be confusing if this code is ever reused for a tool other than mp-excitation-ansatz.
                                    "try using a different value for the option --quantumnumber.";
         throw std::runtime_error(ErrorMessage);
      }
   }

   // Get the null space matrices corresponding to each A-matrix in PsiLeft.
   for (StateComponent C : PsiLinearLeft)
      NullLeftDeque.push_back(NullSpace2(C));
}

std::deque<MatrixOperator>
HEff::operator()(std::deque<MatrixOperator> const& XDeque)
{
   std::deque<StateComponent> BDeque = this->ConstructBDeque(XDeque);

   EF.SetWLower(1, BDeque);

   std::deque<StateComponent> HEffDeque = EF.GetHEff();

   std::deque<StateComponent> TyDeque;
   if (Alpha != 0.0)
   {
      EFTy.SetWLower(1, BDeque);
      TyDeque = EFTy.GetHEff();
   }

   std::deque<MatrixOperator> Result;
   auto NL = NullLeftDeque.begin();
   auto I = HEffDeque.begin();
   auto J = TyDeque.begin();
   while (I != HEffDeque.end())
   {
      Result.push_back(scalar_prod(herm(*NL), *I));
      if (Alpha != 0.0)
      {
         Result.back() += -Alpha * std::conj(ExpIKY) * scalar_prod(herm(*NL), *J);
         ++J;
      }
      ++NL, ++I;
   }

   return Result;
}

std::complex<double>
HEff::Ty(std::deque<MatrixOperator> const& XDeque)
{
   std::deque<StateComponent> BDeque = this->ConstructBDeque(XDeque);

   EFTy.SetWUpper(1, BDeque);
   EFTy.SetWLower(1, BDeque);

   return std::exp(std::complex<double>(0.0, -math_const::pi) * kyFix)
      * inner_prod(EFTy.GetRho({Infinity}, {Infinity}), EFTy.GetElement({Infinity}, {Infinity}).front()[1.0].coefficient(1));
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

std::deque<MatrixOperator>
HEff::InitialGuess() const
{
   std::deque<MatrixOperator> Result;
   auto NL = NullLeftDeque.begin();
   auto CR = PsiLinearRight.begin();
   while (NL != NullLeftDeque.end())
   {
      MatrixOperator C = MakeRandomMatrixOperator((*NL).Basis2(), (*CR).Basis2());
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
      Result.push_back(PackMatrixOperator(MatrixOperator((*NL).Basis2(), (*CR).Basis2())));
      ++NL, ++CR;
   }

   return Result;
}

std::vector<WavefunctionSectionLeft>
HEff::ConstructWindowVec(std::deque<MatrixOperator> XDeque) const
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

   return WindowVec;
}

void
HEff::SetK(double k)
{
   ExpIK = exp(std::complex<double>(0.0, math_const::pi) * k);

   EF.SetExpIKUpper({true}, ExpIK);
   EF.SetExpIKLower({true}, ExpIK);
   if (!StringOp.is_null())
   {
      EFTy.SetExpIKUpper({true}, ExpIK);
      EFTy.SetExpIKLower({true}, ExpIK);
   }
}

void
HEff::SetKY(double ky)
{
   ExpIKY = std::exp(std::complex<double>(0.0, math_const::pi) * (ky + kyFix));
}

PackHEff::PackHEff(HEff* H_)
    : H(H_)
{
   Pack = H->PackInitialize();
   Size = 0;
   for (auto P : Pack)
      Size += P.size();
}

void
PackHEff::operator()(std::complex<double> const* In_, std::complex<double>* Out_) const
{
   std::deque<MatrixOperator> XDeque = this->unpack(In_);
   XDeque = (*H)(XDeque);
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
