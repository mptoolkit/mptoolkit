// -*- C++ -*-
//
// FiniteWavefunction: class to represent a finite-size matrix product wavefunction in canonical form.
//
// The canonical form chosen is Left canonical.  All A-matrices are stored in left-canonical form.
//

#if !defined(MPTOOLKIT_MPS_FINITEWAVEFUNCTION_H)
#define MPTOOLKIT_MPS_FINITEWAVEFUNCTION_H

#include "canonicalwavefunction.h"

class FiniteWavefunction : public CanonicalWavefunctionBase
{
   public:
      FiniteWavefunction();

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
            FiniteMPO const& M, 
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