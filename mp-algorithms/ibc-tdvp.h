// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ibc-tdvp.h
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_IBC_TDVP_H)
#define MPTOOLKIT_MP_ALGORITHMS_IBC_TDVP_H

#include "tdvp.h"
#include "wavefunction/ibc.h"

class IBC_TDVP : public TDVP
{
   public:
      IBC_TDVP() {}

      IBC_TDVP(IBCWavefunction const& Psi_, BasicTriangularMPO const& Ham_,
               std::complex<double> Timestep_, int MaxIter_, double ErrTol_,
               double GMRESTol_, StatesInfo SInfo_, int Verbose_);

      IBCWavefunction Wavefunction() const;

      double GMRESTol;
      InfiniteWavefunctionLeft PsiLeft;
      InfiniteWavefunctionRight PsiRight;
      BasicTriangularMPO HamiltonianLeft;
      BasicTriangularMPO HamiltonianRight;
};

#endif
