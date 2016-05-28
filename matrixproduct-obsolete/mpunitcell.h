// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpunitcell.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
