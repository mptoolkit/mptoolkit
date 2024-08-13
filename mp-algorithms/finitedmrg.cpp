// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp-algorithms/finitedmrg.cpp
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
// $Id$

#include "finitedmrg.h"

FiniteDMRG::FiniteDMRG(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_)
   : DMRG(Verbose_), LastSweepFidelity(0)
{
   this->InitializeLeftOrtho(LinearWavefunction(Psi_.GetSymmetryList(), Psi_.base_begin(), Psi_.base_end()), Ham_, Initial_E(Ham_, Psi_.Basis1()), Initial_F(Ham_, Psi_.Basis2()));

   SweepC = *C;

  this->debug_check_structure();
}

void
FiniteDMRG::AddOrthogonalState(FiniteWavefunctionLeft x)
{
   PANIC("Not implemented");
}

FiniteWavefunctionLeft
FiniteDMRG::Wavefunction() const
{
   return FiniteWavefunctionLeft::Construct(Psi);
}

void
FiniteDMRG::StartSweep()
{
   this->DMRG::StartSweep();

   // Reset the wavfunction for the per-sweep fidelity
   SweepC = *C;
}

void
FiniteDMRG::EndSweep()
{
   LastSweepFidelity = norm_frob(inner_prod(*C, SweepC));
   this->DMRG::EndSweep();
}

void
FiniteDMRG::ModifyLeftBasis(MatrixOperator const& U)
{
   this->DMRG::ModifyLeftBasis(U);
   SweepC = U * SweepC;
}

void
FiniteDMRG::ModifyRightBasis(MatrixOperator const& U)
{
   this->DMRG::ModifyRightBasis(U);
   SweepC = SweepC * U;
}

void
FiniteDMRG::check_structure() const
{
   this->DMRG::check_structure();

   CHECK_EQUAL(C->LocalBasis(), SweepC.LocalBasis());
   CHECK_EQUAL(C->Basis1(), SweepC.Basis1());
   CHECK_EQUAL(C->Basis2(), SweepC.Basis2());
}

void
FiniteDMRG::ShiftLeft(MatrixOperator const& Lambda)
{
   // The projection operator that maps from the SweepC basis to the new basis
   auto L = C; --L;
   SweepC = prod(*L, scalar_prod(SweepC, herm(*C)));
   this->DMRG::ShiftLeft(Lambda);
}

void
FiniteDMRG::ShiftRight(MatrixOperator const& Lambda)
{
   // The projection operator that maps from the SweepC basis to the new basis
   auto R = C; ++R;
   SweepC = prod(scalar_prod(herm(*C), SweepC), *R);
   this->DMRG::ShiftRight(Lambda);
}
