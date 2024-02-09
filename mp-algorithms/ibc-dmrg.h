// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-dmrg.h
//
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_IBC_DMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_IBC_DMRG_H

#include "ibc-tdvp.h"

struct IBC_DMRGSettings : IBC_TDVPSettings
{
};

class IBC_DMRG : public IBC_TDVP
{
   public:
      IBC_DMRG() = default;

      IBC_DMRG(IBCWavefunction const& Psi_, WindowHamiltonian const& Ham_, IBC_DMRGSettings const& Settings_);

      void EvolveCurrentSite(std::complex<double> Tau);
      void IterateLeft(std::complex<double> Tau);
      void IterateRight(std::complex<double> Tau);
};

#endif
