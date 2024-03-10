// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/finitelattice_mpo.h
//
// Copyright (C) 2016-2017 Ian McCulloch <ian@qusim.net>
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

// A FiniteLatticeMPO is an MPO on a FiniteLattice.  It uses a RunLengthCompressed vector
// of GenericMPO, each element being a unit cell.  The elements are stored in triangular format.

#if !defined(MPTOOLKIT_LATTICE_FINITELATTICE_MPO_H)
#define MPTOOLKIT_LATTICE_FINITELATTICE_MPO_H

#include "mpo/generic_mpo.h"
#include "lattice/latticesite.h"

class FiniteLatticeMPO
{
   public:
      typedef BasicFiniteMPO::basis1_type basis1_type;
      typedef BasicFiniteMPO::basis2_type basis2_type;

      FiniteLatticeMPO() = default;

      FiniteLatticeMPO(std::string Description);

      // When we assign a FiniteLatticeMPO, we don't want to copy the description
      // if it has already been set. Hence we need a custom assignment, and therefore
      // other constructors/assignment too.

      FiniteLatticeMPO(FiniteLatticeMPO const&) = default;
      FiniteLatticeMPO(FiniteLatticeMPO&&) = default;

      FiniteLatticeMPO& operator=(FiniteLatticeMPO const& c);
      FiniteLatticeMPO& operator=(FiniteLatticeMPO&& c);

      // returns the total number of sites this operator contains
      int size() const { return Op.size(); }

      // description of the operator
      std::string description() const { return Description; }
      void set_description(std::string s) { Description = std::move(s); }

      // returns true if this is a zero operator
      bool empty() const { return Op.empty(); }
      bool is_null() const { return Op.is_null(); }

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state, and the right basis contains the vacuum.
      // Note that finite MPO's are not irreducible tensor operators as such (they are more like
      // Wigner operators)
      bool is_irreducible() const;

      // returns true if the operator transforms as a rotational invariant, ie
      // it is irreducible in the scalar symmetry sector
      bool is_scalar() const;

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Op.GetSymmetryList(); }

      // If the operator is irreducible, then return the quantum number sector
      QuantumNumbers::QuantumNumber TransformsAs() const { return Op.TransformsAs(); }

      // returns the quantum number in the left basis.  If the right basis is the vacuum
      // then this is also the TransformsAs()
      // precondition: Basis1().size() == 1
      QuantumNumbers::QuantumNumber qn1() const { return Op.qn1(); }

      // returns the quantum number in the right basis.
      // This qill normally be the identity (vacuum)
      QuantumNumbers::QuantumNumber qn2() const { return Op.qn2(); }

      // returns the left-most basis.  This is guaranteed to contain each
      // quantum number at most once.
      basis1_type const& Basis1() const { return Op.Basis1(); }

      // The Basis2() will normally be the identity (vacuum) basis
      basis2_type const& Basis2() const { return Op.Basis2(); }

   private:
      run_length_compressed<GenericMPO> Data;
      BasicFiniteMPO Op;
      std::string Description;

      friend PStream::opstream& operator<<(PStream::opstream& out, FiniteLatticeMPO const& L);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, FiniteLatticeMPO& L);
      friend FiniteLatticeMPO coarse_grain(FiniteLatticeMPO const& Op, int N);
};

std::ostream& operator<<(std::ostream& out, FiniteLatticeMPO const& Op);

FiniteLatticeMPO& operator*=(FiniteLatticeMPO& x, double a);
FiniteLatticeMPO& operator*=(FiniteLatticeMPO& x, std::complex<double> a);

FiniteLatticeMPO& operator+=(FiniteLatticeMPO& x, FiniteLatticeMPO const& y);
FiniteLatticeMPO& operator-=(FiniteLatticeMPO& x, FiniteLatticeMPO const& y);

FiniteLatticeMPO operator+(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);
FiniteLatticeMPO operator-(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

FiniteLatticeMPO operator-(FiniteLatticeMPO const& x);

FiniteLatticeMPO operator*(double a, FiniteLatticeMPO const& x);
FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, double a);
FiniteLatticeMPO operator*(std::complex<double> a, FiniteLatticeMPO const& x);
FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, std::complex<double> a);

FiniteLatticeMPO prod(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y, QuantumNumbers::QuantumNumber const& q);
FiniteLatticeMPO prod(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);
FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
FiniteLatticeMPO dot(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
FiniteLatticeMPO inner(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// cross product (if it exists)
FiniteLatticeMPO cross(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
FiniteLatticeMPO outer(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// power of an operator.  Requires n > 1.
FiniteLatticeMPO pow(FiniteLatticeMPO const& x, int n);

// Exponential operator.
//FiniteLatticeMPO exp(FiniteLatticeMPO const& x);

// Conjugate
FiniteLatticeMPO conj(FiniteLatticeMPO const& x);

// Adjoint
FiniteLatticeMPO adjoint(FiniteLatticeMPO const& x);

// Inverse Adjoint
FiniteLatticeMPO inv_adjoint(FiniteLatticeMPO const& x);

// N into 1 coarse graining
FiniteLatticeMPO coarse_grain(FiniteLatticeMPO const& Op, int N);

// Get initial (1x1) E and F matrices.
StateComponent Initial_E(FiniteLatticeMPO const& m);
StateComponent Initial_F(FiniteMatticeMPO const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(FiniteLatticeMPO const& m, VectorBasis const& B);
StateComponent Initial_F(FiniteLatticeMPO const& m, VectorBasis const& B);

// Optimize the representation
void optimize(FiniteLatticeMPO& Op);

void qr_optimize(FiniteLatticeMPO& Op);

#endif
