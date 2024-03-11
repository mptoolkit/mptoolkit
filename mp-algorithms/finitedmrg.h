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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_FINITEDMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_FINITEDMRG_H

#include "mp-algorithms/dmrg.h"
#include "wavefunction/finitewavefunctionleft.h"

class FiniteDMRG : public DMRG
{
   public:
      FiniteDMRG(FiniteWavefunctionLeft const& Psi_, BasicTriangularMPO const& Ham_, int Verbose_);

      // Adds x to the 'orthogonal set', that we explicitly orthogonalize the
      // wavefunction against
      void AddOrthogonalState(FiniteWavefunctionLeft x);

      // get the current wavefunction
      FiniteWavefunctionLeft Wavefunction() const;

      virtual void check_structure() const;

      //      std::vector<CenterWavefunction> Ortho;              // set of wavefunctions that we want to be
      //      left_right_stack<MatrixOperator> PsiOrthoProjector;  // orthogonal to

};

#endif
