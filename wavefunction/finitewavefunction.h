// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/finitewavefunction.h
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
// FiniteWavefunction: class to represent a finite-size matrix product wavefunction in canonical form.
//
// The canonical form chosen is Left canonical.  All A-matrices are stored in left-canonical form.
//

#if !defined(MPTOOLKIT_MPS_FINITEWAVEFUNCTION_H)
#define MPTOOLKIT_MPS_FINITEWAVEFUNCTION_H

#include "canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_finite_mpo.h"

class AttributeList;

class FiniteWavefunction : public CanonicalWavefunctionBase
{
   public:
      FiniteWavefunction();

      // named constructors

      // construction from a LinearWavefunction (in left-canonical form with lambda
      // matrix on the right)
      static
      FiniteWavefunction ConstructFromOrthogonal(LinearWavefunction const& Psi,
						 MatrixOperator const& Lambda,
						 QuantumNumbers::QuantumNumber const& QShift_,
						 int Verbose = 0);

      // construct and orthogonalize from a LinearWavefunction
      static
      FiniteWavefunctionLeft Construct(LinearWavefunction const& Psi,
				       QuantumNumbers::QuantumNumber const& QShift,
				       int Verbose = 0);

      FiniteWavefunction(FiniteWavefunction const& Psi);

      // construction from a LinearWavefunction
      FiniteWavefunction(LinearWavefunction const& Psi);

      FiniteWavefunction& operator=(FiniteWavefunction const& Psi);

      RealDiagonalOperator lambda_r() const { return this->lambda(this->size()-1); }

      static PStream::VersionTag VersionT;

   private:

};

// calculates the inner product <Psi1|Psi2>
std::complex<double>
overlap(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2);

// calculates the inner product <Psi1|conj(Psi2)>
std::complex<double>
overlap_conj(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2);

// calculates <Psi1|M|Psi2>
std::complex<double>
expectation(FiniteWavefunction const& Psi1,
            BasicFiniteMPO const& M,
            FiniteWavefunction const& Psi2);

double norm_2(FiniteWavefunction const& Psi);

double norm_2_sq(FiniteWavefunction const& Psi);

FiniteWavefunction operator*(double a, FiniteWavefunction const& x);
FiniteWavefunction operator*(std::complex<double> a, FiniteWavefunction const& x);

FiniteWavefunction& operator*=(FiniteWavefunction& psi, double a);
FiniteWavefunction& operator*=(FiniteWavefunction& psi, std::complex<double> a);

FiniteWavefunction conj(FiniteWavefunction const& Psi);

FiniteWavefunction reflect(FiniteWavefunction const& Psi);

HermitianProxy<FiniteWavefunction>
inline
herm(FiniteWavefunction const& x)
{
   return HermitianProxy<FiniteWavefunction>(x);
}


#endif
