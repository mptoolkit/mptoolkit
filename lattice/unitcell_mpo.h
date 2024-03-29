// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/unitcell_mpo.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2023-2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// A UnitCellMPO is a BasicFiniteMPO plus its associated UnitCell.
//
// The UnitCellMPO always starts at a unit cell boundary, and has a size that is
// an integer multiple of the unit cell size.
//
// The Offset parameter of the UnitCellMPO is the site number of the left-hand site of
// the MPO (which is always a multiple of the lattice unit cell size).
//

#if !defined(MPTOOLKIT_LATTICE_UNITCELL_MPO_H)
#define MPTOOLKIT_LATTICE_UNITCELL_MPO_H

#include "mpo/basic_finite_mpo.h"
#include "lattice/latticesite.h"

class UnitCellMPO
{
   public:
      typedef BasicFiniteMPO::basis1_type basis1_type;
      typedef BasicFiniteMPO::basis2_type basis2_type;

      UnitCellMPO() = default;

      UnitCellMPO(SiteListPtrType const& SiteList_, BasicFiniteMPO Op_,
                  LatticeCommute Com_, int Offset_ = 0, std::string Description = "",
		  int CoarseGrain_ = 1);

      // When we assign a UnitCellMPO, we don't want to copy the description
      // if it has already been set. Hence we need a custom assignment, and therefore
      // other constructors/assignment too.

      UnitCellMPO(UnitCellMPO const&) = default;
      UnitCellMPO(UnitCellMPO&&) = default;

      UnitCellMPO& operator=(UnitCellMPO const& c);
      UnitCellMPO& operator=(UnitCellMPO&& c);

      // returns the total number of sites this operator contains
      int size() const { return Op.size(); }

      // returns the size of the unit cell
      int unit_cell_size() const { return SiteList->size(); }

      // if the operator has been coarse-grained, then return the scale factor
      int coarse_grain_factor() const { return CoarseGrain; }

      // The offset - the site number of the first site of the MPO.
      // This must be a multiple of the unit cell size.
      int offset() const { return Offset; }

      // description of the operator
      std::string description() const { return Description; }
      void set_description(std::string s) { Description = std::move(s); }

      // Translates the operator by some number of lattice sites
      void translate(int Shift) { Offset += Shift; }

      // returns true if this is a zero operator
      bool empty() const { return Op.empty(); }
      bool is_null() const { return Op.is_null(); }

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state, and the right basis contains the vacuum.
      // Note that finite MPO's are not irreducible tensor operators as such (they are more like
      // Wigner operators)
      bool is_irreducible() const { return Op.is_irreducible(); }

      // returns true if the operator transforms as a rotational invariant, ie
      // it is irreducible in the scalar symmetry sector
      bool is_scalar() const { return Op.is_scalar(); }

      // returns true if this MPO is the identity operator, that is, a 1x1 MPO that
      // is a product of identity operators.
      //      bool is_identity() const { return Op.is_identity(); }

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Op.GetSymmetryList(); }

      // If the operator is irreducible, then return the quantum number sector
      QuantumNumbers::QuantumNumber TransformsAs() const { return Op.TransformsAs(); }

      // returns the quantum number in the left basis.  If the right basis is the vacuum
      // then this is also the TransformsAs()
      // precondition: Basis1().size() == 1
      QuantumNumbers::QuantumNumber qn1() const { return Op.qn1(); }

      // returns the quantum number in the right basis.
      // This doesn't have to be the vacuum state.
      QuantumNumbers::QuantumNumber qn2() const { return Op.qn2(); }

      // returns the left-most basis.  This is guaranteed to contain each
      // quantum number at most once.
      basis1_type const& Basis1() const { return Op.Basis1(); }

      basis2_type const& Basis2() const { return Op.Basis2(); }

      // Access to the underlying MPO
      BasicFiniteMPO const& MPO() const { return Op; }

      // Access to the underlying MPO
      BasicFiniteMPO& MPO() { return Op; }

      LatticeCommute Commute() const { return Com; }

      SiteListPtrType const& GetSiteList() const { return SiteList; }

      // Extends the MPO with JW-strings on the left and identity on the right
      // so that it covers at least all sites from OtherOffset to OtherOffset+OtherSite.
      void ExtendToCover(int OtherSize, int OtherOffset);

      // Extends the MPO with JW-strings on the left and identity on the right
      // so that it covers complete multiples of OtherSize on both the left and right.
      void ExtendToCoverUnitCell(int OtherSize);

      // returns a representation of the JW string operator as a BasicFiniteMPO
      // acting on a single unit cell
      BasicFiniteMPO GetJWStringUnit() const;

   private:
      SiteListPtrType SiteList;
      BasicFiniteMPO Op;
      LatticeCommute Com;
      int Offset;
      std::string Description;
      int CoarseGrain;

      friend PStream::opstream& operator<<(PStream::opstream& out, UnitCellMPO const& L);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, UnitCellMPO& L);
      friend UnitCellMPO coarse_grain(UnitCellMPO const& Op, int N);
};

std::ostream& operator<<(std::ostream& out, UnitCellMPO const& Op);

UnitCellMPO& operator*=(UnitCellMPO& x, double a);
UnitCellMPO& operator*=(UnitCellMPO& x, std::complex<double> a);

UnitCellMPO& operator+=(UnitCellMPO& x, UnitCellMPO const& y);
UnitCellMPO& operator-=(UnitCellMPO& x, UnitCellMPO const& y);

UnitCellMPO operator+(UnitCellMPO const& x, UnitCellMPO const& y);
UnitCellMPO operator-(UnitCellMPO const& x, UnitCellMPO const& y);

UnitCellMPO operator-(UnitCellMPO const& x);

UnitCellMPO operator*(double a, UnitCellMPO const& x);
UnitCellMPO operator*(UnitCellMPO const& x, double a);
UnitCellMPO operator*(std::complex<double> a, UnitCellMPO const& x);
UnitCellMPO operator*(UnitCellMPO const& x, std::complex<double> a);

UnitCellMPO prod(UnitCellMPO const& x, UnitCellMPO const& y, QuantumNumbers::QuantumNumber const& q);
UnitCellMPO prod(UnitCellMPO const& x, UnitCellMPO const& y);
UnitCellMPO operator*(UnitCellMPO const& x, UnitCellMPO const& y);

UnitCellMPO commutator(UnitCellMPO const& x, UnitCellMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
UnitCellMPO dot(UnitCellMPO const& x, UnitCellMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
UnitCellMPO inner(UnitCellMPO const& x, UnitCellMPO const& y);

// cross product (if it exists)
UnitCellMPO cross(UnitCellMPO const& x, UnitCellMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
UnitCellMPO outer(UnitCellMPO const& x, UnitCellMPO const& y);

// project a (reducible) operator onto an irreducible component
UnitCellMPO project(UnitCellMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.
UnitCellMPO pow(UnitCellMPO const& x, int n);

// Exponential operator.
UnitCellMPO exp(UnitCellMPO const& x);

// Absolute value of an operator.
UnitCellMPO abs(UnitCellMPO const& x);

// Conjugate
UnitCellMPO conj(UnitCellMPO const& x);

// Adjoint
UnitCellMPO adjoint(UnitCellMPO const& x);

// Inverse Adjoint
UnitCellMPO inv_adjoint(UnitCellMPO const& x);

UnitCellMPO gauge_flip(UnitCellMPO const& x);

// translate - shift a UnitCellMPO by some number of sites.
// This can be positive or negative but MUST be a multiple of
// the unit cell size.  TODO: relax this restriction as long as the
// SiteList is invariant under the shift
UnitCellMPO translate(UnitCellMPO x, int Sites);

// N into 1 coarse graining
UnitCellMPO coarse_grain(UnitCellMPO const& Op, int N);

// Constructs an identity MPO from a given unit cell
UnitCellMPO MakeIdentityFrom(UnitCellMPO const& x);

// Free versions of ExtendToCover and ExtendToCoverUnitCell
inline
UnitCellMPO ExtendToCover(UnitCellMPO const& Op, int OtherSize, int OtherOffset)
{
   UnitCellMPO Result(Op);
   Result.ExtendToCover(OtherSize, OtherOffset);
   return Result;
}

inline
UnitCellMPO ExtendToCoverUnitCell(UnitCellMPO const& Op, int OtherSize)
{
   UnitCellMPO Result(Op);
   Result.ExtendToCoverUnitCell(OtherSize);
   return Result;
}

inline
UnitCellMPO coarse_grain(UnitCellMPO const& Op, int N)
{
   return UnitCellMPO(Op.SiteList, coarse_grain(ExtendToCoverUnitCell(Op,N).Op, N), Op.Com,
                      numerics::divp(Op.Offset, N).quot, Op.Description,
		      Op.CoarseGrain*N);
}

// Optimize the representation - in this case we simply forward to the BasicFiniteMPO representation
void optimize(UnitCellMPO& Op);

void qr_optimize(UnitCellMPO& Op);

#endif
