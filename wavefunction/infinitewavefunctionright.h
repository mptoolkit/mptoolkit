// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/infinitewavefunctionright.h
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
// InfiniteWavefunctionRight: class to represent a linear `infinite' matrix product wavefunction
// in right canonical form.
//

#if !defined(MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONRIGHT_H)
#define MPTOOLKIT_WAVEFUNCTION_INFINITEWAVEFUNCTIONRIGHT_H

#include "wavefunction/canonicalwavefunction.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/basic_finite_mpo.h"
#include "mpo/product_mpo.h"
#include <tuple>

class InfiniteWavefunctionLeft;

// class to represent an infinite wavefunction in the left canonical basis
class InfiniteWavefunctionRight : public CanonicalWavefunctionBase
{
   public:
      InfiniteWavefunctionRight() = default;
      InfiniteWavefunctionRight(InfiniteWavefunctionRight const& x) = default;
      InfiniteWavefunctionRight(InfiniteWavefunctionRight&& x) = default;
      InfiniteWavefunctionRight& operator=(InfiniteWavefunctionRight const& Psi) = default;
      InfiniteWavefunctionRight& operator=(InfiniteWavefunctionRight&& Psi) = default;

      // Construction from an InfiniteWavefunctionLeft
      explicit InfiniteWavefunctionRight(InfiniteWavefunctionLeft const& Psi);

      QuantumNumber qshift() const { return QShift; }

      // Rotates the wavefunction to the left, by taking the left-most site and moving
      // it to the right
      void rotate_left(int Count);

      // Rotates the wavefunction to the left, by taking the right-most site and moving
      // it to the left
      void rotate_right(int Count);

      void SetDefaultAttributes(AttributeList& A) const;

      static std::string Type;

      static PStream::VersionTag VersionT;

      friend PStream::ipstream& operator>>(PStream::ipstream& in,
                                           InfiniteWavefunctionRight& Psi);
      friend PStream::opstream& operator<<(PStream::opstream& out,
                                           InfiniteWavefunctionRight const& Psi);
      friend void read_version(PStream::ipstream& in, InfiniteWavefunctionRight& Psi, int Version);

      void check_structure() const;
      void debug_check_structure() const;

   private:
      void Initialize(MatrixOperator const& Lambda, LinearWavefunction const& Psi);

      QuantumNumber QShift;

      // All functions that can modify the internal representation but preserve the canonical form
      // are friend functions.  This is so that we have a central list of such functions,
      // so can update them if the class changes.
      friend void inplace_reflect(InfiniteWavefunctionRight& Psi);
      friend void inplace_conj(InfiniteWavefunctionRight& Psi);
      friend void inplace_qshift(InfiniteWavefunctionRight& Psi, QuantumNumbers::QuantumNumber const& Shift);
      friend InfiniteWavefunctionRight wigner_project(InfiniteWavefunctionRight const& Psi,
                                                      SymmetryList const& FinalSL);
      friend InfiniteWavefunctionRight ReorderSymmetry(InfiniteWavefunctionRight const& Psi,
                                                       SymmetryList const& NewSL);
      friend InfiniteWavefunctionRight repeat(InfiniteWavefunctionRight const& Psi, int Count);
};

class InfiniteWavefunctionLeft;

// Convert a, infinite wavefunction to right-canonical form,
// and returns the Lambda matrix on the right-hand-side.
RealDiagonalOperator
right_canonicalize(LinearWavefunction& Psi, QuantumNumbers::QuantumNumber const& QShift);

// Extend the unit cell of the wavefunction by repeating it Count number of times
InfiniteWavefunctionRight repeat(InfiniteWavefunctionRight const& Psi, int Count);

// returns a linear wavefunction that is in pure right-orthogonal form.
// This is a very fast operation that only manipulates pvalue_handle objects.
std::pair<RealDiagonalOperator, LinearWavefunction>
get_right_canonical(InfiniteWavefunctionRight const& Psi);

// returns a linear wavefunction that is in pure left-orthogonal form.
// The asymmetry in the return type between get_right_canonical and get_left_canonical
// is because for a right-canonical wavefunction, the Lambda matrix is already in diagonal form,
// whereas converting a right-canonical wavefunction into a left-canonical wavefunction gives
// an additional unitary matrix, that we cannot wrap around to the other end of the wavefunction
// as it would change the basis there.
// The wavefunction isn't quite useable as-is, the additional unitary needs to be incorporated.
// There are two ways to do this:
// Keeping the basis at the unit cell boundary fixed, transform the lambda matrix as
// herm(U)*lambda*(U) and transform Psi as Psi*U.  This makes lambda no longer a diagonal operator.
// The other alternative is to keep lambda as diagonal, but change the basis at the unit cell
// boundary.  To do this, transform Psi as U*Psi.
std::tuple<LinearWavefunction, RealDiagonalOperator, MatrixOperator>
get_left_canonical(InfiniteWavefunctionRight const& Psi);

// function to extract the local basis (as a vector of BasisList) from a wavefunction
// TODO: this implementation isn't terrible, but isn't so elegant either!
std::vector<BasisList>
inline
ExtractLocalBasis(InfiniteWavefunctionRight const& Psi)
{
   return ExtractLocalBasis(get_right_canonical(Psi).second);
}

// calculates the overlap of two iMPS, per unit cell.
// The eigenvector can be in any allowable symmetry sector.
// x and y must have the same size
std::pair<std::complex<double>, StateComponent>
overlap(InfiniteWavefunctionRight const& x,  InfiniteWavefunctionRight const& y,
        QuantumNumbers::QuantumNumber const& Sector,
        int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// This version allows the wavefunctions and operator to have different sizes.
// The overlap is returned as a quantity per length, which is the lowest
// common multiple of x.size(), y.size(), StringOp.size()
// The length is returned as the second component of the tuple
std::tuple<std::complex<double>, int, StateComponent>
overlap(InfiniteWavefunctionRight const& x, ProductMPO const& StringOp,
        InfiniteWavefunctionRight const& y,
        QuantumNumbers::QuantumNumber const& Sector,
        int Iter = 20, double Tol = 1E-12, int Verbose = 0);

// Reflect a wavefunction in place
void inplace_reflect(InfiniteWavefunctionRight& Psi);

// Conjugate a wavefunction in place
void inplace_conj(InfiniteWavefunctionRight& Psi);

// Spatial reflection of a wavefunction
InfiniteWavefunctionLeft reflect(InfiniteWavefunctionRight const& Psi);

// version of reflect where we apply a local operator also
//InfiniteWavefunctionRight reflect(InfiniteWavefunctionLeft const& Psi,
// std::vector<SimpleOperator> const& Op);

inline
void
InfiniteWavefunctionRight::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

#endif
