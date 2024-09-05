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

      // get the current wavefunction. This could be called at any time.
      FiniteWavefunctionLeft Wavefunction() const;

      void StartSweep();
      void EndSweep();

   protected:
      virtual void ShiftLeft(MatrixOperator const& Lambda);
      virtual void ShiftRight(MatrixOperator const& Lambda);

      virtual void ModifyLeftBasis(MatrixOperator const& U);
      virtual void ModifyRightBasis(MatrixOperator const& U);

   public:
      virtual void check_structure() const;

      // Statistics which are valid at the end of the sweep
      double LastSweepFidelity;      // Overlap of the wavefunction before and after the previous sweep


      // The wavefunction from the start of the sweep, projected into the basis of C
      StateComponent SweepC;

      //      std::vector<CenterWavefunction> Ortho;              // set of wavefunctions that we want to be
      //      left_right_stack<MatrixOperator> PsiOrthoProjector;  // orthogonal to

};

#endif
