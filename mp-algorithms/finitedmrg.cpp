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
   : DMRG(Verbose_)
{
   this->InitializeLeftOrtho(LinearWavefunction(Psi_.GetSymmetryList(), Psi_.base_begin(), Psi_.base_end()), Ham_, Initial_E(Ham_, Psi_.Basis1()), Initial_F(Ham_, Psi_.Basis2()));

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
FiniteDMRG::check_structure() const
{
   this->DMRG::check_structure();
}
