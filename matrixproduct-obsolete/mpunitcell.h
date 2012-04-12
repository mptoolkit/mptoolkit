// -*- C++ -*- $id$

// defines MPUnitCell, for the unit cell of an infinite system.
// This refines a LinearWavefunction such that the Basis1() is always
// the same as the Basis2(), and incorporates a quantum number shift 
// for Abelian quantum numbers.

#if !defined(MPUNITCELL_H_8489HUH5789TY789HPO98)
#define MPUNITCELL_H_8489HUH5789TY789HPO98

#include "linearwavefunction.h"

class MPUnitCell
{
   public:
      typedef MPStateComponent value_type;
   
      MPUnitCell() {}

      MPUnitCell(LinearWavefunction const& Psi, QuantumNumber const& QShift);

      // Returns the basis that represents the wrap-around of the unit cell
      VectorBasis Basis() const;

      // for compatibility with other MPS types.  Returns Basis()
      VectorBasis Basis1() const;
      // for compatibility with other MPS types.  Returns Basis()
      VectorBasis Basis2() const;

   private:
      LinearWavefunction Psi_;
      QuantumNumber QShift_;
};

HermitianProxy<MPUnitCell>
inline
herm(MPUnitCell const& x)
{
   return HermitianProxy<MPUnitCell>(x);
}

// from left to right, calculates the action of the transfer operator R = A^\dagger m B
MatrixOperator
operator_prod(HermitianProxy<MPUnitCell> const& A, 
              MatrixOperator const& m, 
              MPUnitCell const& B);

// from right to left, calculates the action of the transfer operator R = A m B^\dagger
MatrixOperator
operator_prod(MPUnitCell const& A,
              MatrixOperator const& m, 
              HermitianProxy<MPUnitCell> const& B);





#endif
