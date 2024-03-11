// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/finite_lattice_mpo.h
//
// Copyright (C) 2014-2017 Ian McCulloch <ian@qusim.net>
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
// a FiniteLatticeMPO is a BasicFiniteMPO that also knows about the
// unit cell, and it also knows what its commutation Jordan-Wigner string looks like.
// This means that it is possible to extend a FiniteLatticeMPO by adding a unit
// cell to either the left or the right.

#if !defined(MPTOOLKIT_MPO_FINITE_LATTICE_MPO_H)
#define MPTOOLKIT_MPO_FINITE_LATTICE_MPO_H

#include "finite_mpo.h"
#include "siteoperator/unitcell.h"

class FiniteLatticeMPO
{
   public:
      FiniteLatticeMPO();

      FiniteLatticeMPO(FiniteLatticeMPO const& Other);

      FiniteLatticeMPO& operator=(FiniteLatticeMPO const& Other);

      // Construction via components
      FiniteLatticeMPO(UnitCell const& UC, BasicFiniteMPO const& JW, BasicFiniteMPO const& Op);

      // returns the number of sites that this operator is defined over.
      // This is always a multiple of the unit cell size.
      int size() const;

      // returns the Jordan-Wigner string associated with this operator,
      // This is always a 1x1 MPO, with the same size as the unit cell.
      BasicFiniteMPO const& JWString() const { return JWString_; }

      // returns the Jordan-Wigner string associated with this operator,
      // extended to be Size sites (which must be a multiple of the UnitCellSize)
      BasicFiniteMPO JWString(int Size) const;

      // returns the UnitCell of this operator
      UnitCell const& GetUnitCell() const { return UnitCell_; }

      // shorthand for GetUnitCell().size()
      int UnitCellSize() const { return UnitCell_.size(); }

      // implicit conversion to BasicFiniteMPO
      operator BasicFiniteMPO const&() const { return Operator_; }

      // named conversion to BasicFiniteMPO
      BasicFiniteMPO const& AsBasicFiniteMPO() const { return Operator_; }

      // Multiplies this operator by the Jordan-Wigner string of f, and returns
      // the BasicFiniteMPO.  This corresponds to the operator obtained on the left hand
      // side when multiplying by operator f that acts on a unit cell off the right hand
      // side of this operator, ie *this \oprod f
      BasicFiniteMPO ApplyJW(FiniteLatticeMPO const& f);

      // Extends the operator to be Size sites (must be a multiple of the UnitCellSize)
      // by adding identity operators on the right hand side.
      BasicFiniteMPO AsBasicFiniteMPO(int Size) const;

   private:
      UnitCell UnitCell_;
      BasicFiniteMPO JWString_;
      BasicFiniteMPO Operator_;

      friend PStream::opstream&
         operator<<(PStream::opstream& out, FiniteLatticeMPO const& op);
      friend PStream::ipstream&
         operator>>(PStream::ipstream& in, FiniteLatticeMPO& op);

      friend FiniteLatticeMPO& operator*=(FiniteLatticeMPO& x, double a);
      friend FiniteLatticeMPO& operator*=(FiniteLatticeMPO& x, std::complex<double> a);

      friend FiniteLatticeMPO& operator+=(FiniteLatticeMPO& x, FiniteLatticeMPO const& y);
      friend FiniteLatticeMPO& operator-=(FiniteLatticeMPO& x, FiniteLatticeMPO const& y);
};

// Returns the identity operator for the given unit cell
BasicFiniteMPO
identity_mpo(UnitCell const& c);

// Returns the identity operator for the given unit cell, extended
// to Size sites, which must be a multiple of the unit cell size.
BasicFiniteMPO
identity_mpo(UnitCell const& c, int Size);

// extends a FiniteLatticeMPO to Size sites, by adding identity operators
// on the right hand side.  Size must be a multiple of f.UnitCellSize().
FiniteLatticeMPO
extend_to(FiniteLatticeMPO const& f, int Size);

// extends a FiniteLatticeMPO by adding LeftSites sites on the left (as the J-W string),
// and RightSites sites to the right (as identity operators).
FiniteLatticeMPO
extend(int LeftSites, FiniteLatticeMPO const& f, int RightSites);

// operations

PStream::opstream& operator<<(PStream::opstream& out, FiniteLatticeMPO const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, FiniteLatticeMPO& op);

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

FiniteLatticeMPO prod(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y,
                      QuantumNumbers::QuantumNumber const& q);
FiniteLatticeMPO prod(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);
FiniteLatticeMPO operator*(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
FiniteLatticeMPO dot(FiniteLatticeMPO const& x, FiniteLatticeMPO const& y);

// project a (reducible) quantum number onto an irreducible component
FiniteLatticeMPO project(FiniteLatticeMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.  Only useful for n small!
FiniteLatticeMPO pow(FiniteLatticeMPO const& x, int n);

// Conjugate
FiniteLatticeMPO conj(FiniteLatticeMPO const& x);

// Adjoint
FiniteLatticeMPO adjoint(FiniteLatticeMPO const& x);

#endif
