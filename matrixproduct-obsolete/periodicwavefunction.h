// -*- C++ -*- $Id$
//
// PeriodicWavefunction: main class to represent a 
// translationally invariant periodic matrix product wavefunction on a finite lattice.
// Currently, the unit cell is a single site.
//

#if !defined(PERIODICWAVEFUNCTION_H_FUIYT49786Y709)
#define PERIODICWAVEFUNCTION_H_FUIYT49786Y709

#include "mpstate.h"
#include "pheap/pvalueptr.h"
#include "pheap/pvalueiterator.h"
#include "interface/attributes.h"
#include "density.h"
#include "linearoperator.h"

class PeriodicWavefunction
{
   public:
      typedef MPStateComponent value_type;
      typedef value_type::value_type operator_type;  // typically this will be MatrixOperator

      PeriodicWavefunction() {}

      // Construct a periodic wavefunction from a given MPStateComponent.  Assumed to transform as a scalar.
      PeriodicWavefunction(int Size, MPStateComponent const& Data);

      PeriodicWavefunction(int Size, MatrixOperator const& Q, MPStateComponent const& Data);

      PeriodicWavefunction(int Size, BasisList const& LocalBasis, VectorBasis const& Basis);

      PeriodicWavefunction(int Size, BasisList const& LocalBasis, VectorBasis const& Basis,
			   QuantumNumber const& Trans);

      QuantumNumber TransformsAs() const { return adjoint(Q_.TransformsAs()); }

      SymmetryList GetSymmetryList() const { return Data_.GetSymmetryList(); }

      std::size_t size() const { return Size_; }
      void set_size(int Size) { Size_ = Size; }

      BasisList SiteBasis() const { return Data_.SiteBasis(); }
      VectorBasis Basis() const { return Data_.Basis1(); }
 
      MatrixOperator const& Q() const { return Q_; }
      MatrixOperator & Q() { return Q_; }

      MPStateComponent const& Data() const { return Data_; }
      MPStateComponent& Data() { return Data_; }

      // Direct access to the matrix components of a given local basis state
      operator_type const& operator[](int s) const { return Data_[s]; }
      operator_type& operator[](int s) { return Data_[s]; }

   private:
      int Size_;
      MatrixOperator Q_;
      MPStateComponent Data_;

   friend PStream::opstream& operator<<(PStream::opstream& out, PeriodicWavefunction const& psi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, PeriodicWavefunction& psi);
};

std::ostream& operator<<(std::ostream& out, PeriodicWavefunction const& psi);

#endif
