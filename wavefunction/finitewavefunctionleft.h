// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/finitewavefunctionleft.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
//
// FiniteWavefunctionLeft: class to represent a finite-size matrix product wavefunction
// in left canonical form.
//

#if !defined(MPTOOLKIT_MPS_FINITEWAVEFUNCTIONLEFT_H)
#define MPTOOLKIT_MPS_FINITEWAVEFUNCTIONLEFT_H

#include "wavefunction/canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_finite_mpo.h"

class AttributeList;

class FiniteWavefunctionLeft : public CanonicalWavefunctionBase
{
   public:
      // Compiler generated constructors are OK

      // named constructors

      // construction from a LinearWavefunction that is right-orthogonalized
      // Psi must be normalized, but is multiplied by an amplitude factor a
      static
      FiniteWavefunctionLeft ConstructFromRightOrthogonal(LinearWavefunction Psi,
							  std::complex<double> a,
							  int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction
      static
      FiniteWavefunctionLeft Construct(LinearWavefunction Psi,
				       int Verbose = 0);

      void SetDefaultAttributes(AttributeList& A) const;

      static PStream::VersionTag VersionT;

      friend FiniteWavefunctionLeft& operator*=(FiniteWavefunctionLeft& psi, double a);
      friend FiniteWavefunctionLeft& operator*=(FiniteWavefunctionLeft& psi,
						std::complex<double> a);
      friend void inplace_conj(FiniteWavefunctionLeft& Psi);
      friend void inplace_reflect(FiniteWavefunctionLeft& Psi);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, FiniteWavefunctionLeft& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out,
                                           FiniteWavefunctionLeft const& Psi);

      friend FiniteWavefunctionLeft
      wigner_project(FiniteWavefunctionLeft const& Psi, SymmetryList const& FinalSL);
      friend FiniteWavefunctionLeft ReorderSymmetry(FiniteWavefunctionLeft const& Psi, SymmetryList const& NewSL);

      static std::string const Type;
};

// calculates the inner product <Psi1|Psi2>
std::complex<double>
overlap(FiniteWavefunctionLeft const& Psi1, FiniteWavefunctionLeft const& Psi2);

// calculates the inner product <Psi1|conj(Psi2)>
std::complex<double>
overlap_conj(FiniteWavefunctionLeft const& Psi1, FiniteWavefunctionLeft const& Psi2);

// calculates <Psi1|M|Psi2>
std::complex<double>
expectation(FiniteWavefunctionLeft const& Psi1,
            BasicFiniteMPO const& M,
            FiniteWavefunctionLeft const& Psi2);

double norm_2(FiniteWavefunctionLeft const& Psi);

double norm_2_sq(FiniteWavefunctionLeft const& Psi);

inline
void normalize(FiniteWavefunctionLeft& Psi)
{
   Psi *= 1.0 / norm_2(Psi);
}

FiniteWavefunctionLeft operator*(double a, FiniteWavefunctionLeft x);
FiniteWavefunctionLeft operator*(FiniteWavefunctionLeft x, double a);
FiniteWavefunctionLeft operator*(std::complex<double> a, FiniteWavefunctionLeft x);
FiniteWavefunctionLeft operator*(FiniteWavefunctionLeft x, std::complex<double> a);

FiniteWavefunctionLeft& operator*=(FiniteWavefunctionLeft& psi, double a);
FiniteWavefunctionLeft& operator*=(FiniteWavefunctionLeft& psi, std::complex<double> a);

FiniteWavefunctionLeft conj(FiniteWavefunctionLeft const& Psi);

FiniteWavefunctionLeft reflect(FiniteWavefunctionLeft const& Psi);

void inplace_reflect(FiniteWavefunctionLeft& Psi);

void inplace_conj(FiniteWavefunctionLeft& Psi);

#endif
