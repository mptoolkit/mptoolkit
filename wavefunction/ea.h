// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/ea.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// Excitation anzatz wavefunction.
//
// This is essentially the same as an IBC wavefunction but with three
// differences:
//
// 1. There are multiple windows, corresponding to each site in the unit cell
// (so we assume that Left, WindowVec and Right have the same size, which could
// be checked or relaxed: TODO). This is more efficient than writing a big
// window which combines all of the components.
//
// 2. The EA wavefunction is essentially equal to the momentum superposition of
// an IBC with momentum k:
// EA = \sum_n exp(i*k*n) T^n IBC,
// where T is the unit-cell translation operator. We store the value of the
// phase exp(i*k) as ExpIK.
//
// 3. By construction, EA wavefunctions are orthogonal to the ground state,
// which is enforced explicitly by writing the left site of the window as the
// null space of the ground state A-matrix. Here, we assume that the windows
// are orthogonal to the ground state which could be checked explictly in the
// future (TODO). In the case that we may want to have an EA wavefunction with
// a component in the direction of the ground state, we still keep the windows
// orthogonal, but track this component with the GSOverlap term.

#if !defined(MPTOOLKIT_WAVFUNCTION_EA_H)
#define MPTOOLKIT_WAVFUNCTION_EA_H

#include "infinitewavefunctionleft.h"
#include "infinitewavefunctionright.h"
#include "canonicalwavefunction.h"
#include "ibc.h"
#include "lattice/unitcell_mpo.h"

class EAWavefunction
{
   public:
      EAWavefunction();

      EAWavefunction(EAWavefunction const& Psi) = default;

      EAWavefunction(InfiniteWavefunctionLeft const& Left_,
                     std::vector<WavefunctionSectionLeft> const& WindowVec_,
                     InfiniteWavefunctionRight const& Right_,
                     std::complex<double> ExpIK = 1.0,
                     std::complex<double> GSOverlap = 0.0);

      SymmetryList GetSymmetryList() const { return WindowVec.front().GetSymmetryList(); }

      // Assume that each window has the same size.
      int window_size() const { return WindowVec.front().size(); }

      void SetDefaultAttributes(AttributeList& A) const;

      static std::string const Type;

      static PStream::VersionTag VersionT;

      void check_structure() const;
      void debug_check_structure() const;

      // private:

      InfiniteWavefunctionLeft Left;
      std::vector<WavefunctionSectionLeft> WindowVec;
      InfiniteWavefunctionRight Right;

      std::complex<double> ExpIK;
      std::complex<double> GSOverlap;

      friend void inplace_reflect(EAWavefunction& Psi);
      friend void inplace_conj(EAWavefunction& Psi);

      friend EAWavefunction wigner_project(EAWavefunction const& Psi,
                                           SymmetryList const& FinalSL);
      friend EAWavefunction ReorderSymmetry(EAWavefunction const& Psi,
                                            SymmetryList const& NewSL);

      friend PStream::opstream& operator<<(PStream::opstream& out, EAWavefunction const& Psi);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, EAWavefunction& Psi);

};

// Reflect a wavefunction in place
void inplace_reflect(EAWavefunction& Psi);

// Conjugate a wavefunction in place
void inplace_conj(EAWavefunction& Psi);

// Spatial reflection of a wavefunction
EAWavefunction reflect(EAWavefunction const& Psi);

#endif
