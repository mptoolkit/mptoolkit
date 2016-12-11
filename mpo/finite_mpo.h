// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/finite_mpo.h
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
//
// Matrix product operator defined on finite support.
// The boundary states are normally one dimensional
// (we also allow the operator to be reducible,
// representing a sum of quantum number components, in which case the Basis1() will have dimension > 1).
// We used to require that the Basis2() was a scalar, but no longer,
// this isn't possible for extracted components of triangular or
// generic operators.  But we can always do a delta_shift to give a scalar (can we??!?)

#if !defined(FINITE_MPO_H_JDCHJKEHY589758YUER89H489)
#define FINITE_MPO_H_JDCHJKEHY589758YUER89H489

#include "generic_mpo.h"
#include "lattice/latticesite.h"

class FiniteMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type;
      typedef data_type::const_iterator const_iterator;
      typedef data_type::iterator       iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      FiniteMPO() {}

      FiniteMPO(FiniteMPO const& Other) : Data(Other.Data) {}

      // Removed this constructor because it doesn't make much sense to define a FiniteMPO
      // without specifying the LatticeCommute
      //      explicit FiniteMPO(int Size) : Data(Size) {}

      explicit FiniteMPO(int Size) : Data(Size) {}

      // Construction from a generic MPO.  The generic MPO must already be in finite form.
      explicit FiniteMPO(GenericMPO const& Other);

      FiniteMPO& operator=(FiniteMPO const& Other) { Data = Other.Data; return *this; }

      // returns the total number of sites this operator contains
      int size() const { return Data.size(); }

      // returns true if this is a zero operator
      bool empty() const { return Data.empty() || Data.front().Basis1().size() == 0; }
      bool is_null() const { return Data.empty() || Data.front().Basis1().size() == 0; }

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state, and the right basis contains the vacuum.
      // Note that finite MPO's are not irreducible tensor operators as such (they are more like
      // Wigner operators)
      bool is_irreducible() const;

      // returns true if the operator transforms as a rotational invariant, ie
      // it is irreducible in the scalar symmetry sector
      bool is_scalar() const;

      // precondition: is_irreducible
      // WARNING: Use qn1() instead of TransformsAs() where appropriate.
      // For debugging purposes, to audit usage of TransformsAs, we have an include guard
#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
      QuantumNumbers::QuantumNumber TransformsAs() const;
#endif

      // returns the quantum number in the left basis.  If the right basis is the vacuum
      // then this is also the TransformsAs()
      // precondition: Basis1().size() == 1
      QuantumNumbers::QuantumNumber qn1() const;

      // returns the quantum number in the right basis.
      // This doesn't have to be the vacuum state.
      QuantumNumbers::QuantumNumber qn2() const;

      // returns the left-most basis.  This is guaranteed to contain each
      // quantum number at most once.
      basis1_type const& Basis1() const { return Data.front().Basis1(); }

      basis2_type const& Basis2() const { return Data.back().Basis2(); }

      value_type& operator[](int n) { return Data[n]; }
      value_type const& operator[](int n) const { return Data[n]; }

      // iterate over the MPOpComponents at each site
      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      value_type& front() { return Data.front(); }
      value_type const& front() const { return Data.front(); }

      value_type& back() { return Data.back(); }
      value_type const& back() const { return Data.back(); }

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Data.front().GetSymmetryList(); }

      // implicit conversion to a const GenericMPO
      operator GenericMPO const&() const { return Data; }

      // Return the local basis at the n'th site
      BasisList const& LocalBasis1(int n) const
      { return Data.LocalBasis1(n); }
      BasisList const& LocalBasis2(int n) const
      { return Data.LocalBasis2(n); }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const
      { return Data.LocalBasis1List(); }
      std::vector<BasisList> LocalBasis2List() const
      { return Data.LocalBasis2List(); }

      // direct access to the GenericMPO
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

      void check_structure() const { Data.check_structure(); }
      void debug_check_structure() const { Data.debug_check_structure(); }

      // Make an identity MPO over the given unit cell basis
      static FiniteMPO make_identity(std::vector<BasisList> const& Basis);
      static FiniteMPO make_identity(std::vector<BasisList> const& Basis,
                                     QuantumNumber const& q);

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, FiniteMPO const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, FiniteMPO& op);

// Returns the MPO that is Op1 \otimes Op2.
// PRECONDITION: Op1.Basis2() == Op2.Basis1()
FiniteMPO join(FiniteMPO const& Op1, FiniteMPO const& Op2);

// Repeats Op Count number of times, Op^{\oprod Count}.
// PRECONDITION: Op.Basis2() == Op.Basis1()
FiniteMPO repeat(FiniteMPO const& Op, int Count);

FiniteMPO& operator*=(FiniteMPO& x, double a);
FiniteMPO& operator*=(FiniteMPO& x, std::complex<double> a);

FiniteMPO& operator+=(FiniteMPO& x, FiniteMPO const& y);
FiniteMPO& operator-=(FiniteMPO& x, FiniteMPO const& y);

FiniteMPO operator+(FiniteMPO const& x, FiniteMPO const& y);
FiniteMPO operator-(FiniteMPO const& x, FiniteMPO const& y);

FiniteMPO operator-(FiniteMPO const& x);

FiniteMPO operator*(double a, FiniteMPO const& x);
FiniteMPO operator*(FiniteMPO const& x, double a);
FiniteMPO operator*(std::complex<double> a, FiniteMPO const& x);
FiniteMPO operator*(FiniteMPO const& x, std::complex<double> a);

FiniteMPO prod(FiniteMPO const& x, FiniteMPO const& y, QuantumNumbers::QuantumNumber const& q);
FiniteMPO prod(FiniteMPO const& x, FiniteMPO const& y);
FiniteMPO operator*(FiniteMPO const& x, FiniteMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
FiniteMPO dot(FiniteMPO const& x, FiniteMPO const& y);

// cross product (if it exists)
FiniteMPO cross(FiniteMPO const& x, FiniteMPO const& y);

// Helper function for the coupling coefficient for the outer() function
double outer_coefficient(int degree_x, int degree_y, int degree_q);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
FiniteMPO outer(FiniteMPO const& x, FiniteMPO const& y);

// project a (reducible) operator onto an irreducible component
FiniteMPO project(FiniteMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.
FiniteMPO pow(FiniteMPO const& x, int n);

// Exponential operator.
FiniteMPO exp(FiniteMPO const& x);

// Conjugate
FiniteMPO conj(FiniteMPO const& x);

// Adjoint
FiniteMPO adjoint(FiniteMPO const& x);

// Inverse Adjoint
FiniteMPO inv_adjoint(FiniteMPO const& x);

// optimize the representation
void optimize(FiniteMPO& Op);

// optimize the representation using QR decomposition
void qr_optimize(FiniteMPO& Op);

// completely coarse-grain the MPO into a simple operator.
// The dimensions of this operator are exponentially big in the number of sites
// in x, so be careful!
// For non-abelian symmetries, this coarse-grain occurs from left to right.
SimpleRedOperator coarse_grain(FiniteMPO const& x);

// Does a N-1 coarse graining of an operator.  The length must be a multiple of N
FiniteMPO coarse_grain(FiniteMPO const& Op, int N);

// The opposite of coarse_grain - decompose an operator acting on the entire Hilbert space
// into a FiniteMPO
FiniteMPO fine_grain(SimpleOperator const& x,
                     std::vector<BasisList> const& LocalBasis1,
                     std::vector<BasisList> const& LocalBasis2);

// Make an identity operator that acts on the same local Hilbert space as x
FiniteMPO
MakeIdentityFrom(FiniteMPO const& x);

// Make an identity operator that acts on the same local Hilbert space as x,
// with the given quantum number in the auxiliary basis
FiniteMPO
MakeIdentityFrom(FiniteMPO const& x, QuantumNumber const& q);

// output to a stream
std::ostream& operator<<(std::ostream& out, FiniteMPO const& x);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(FiniteMPO const& Op, std::ostream& out, double UnityEpsilon);

inline
void print_structure(FiniteMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

// returns the FiniteMPO for the identity operator acting on the unit cell
// with quantum number q in the auxiliary basis
FiniteMPO identity_mpo(SiteListType const& SiteList, QuantumNumbers::QuantumNumber const& q);

// returns the FiniteMPO for the identity operator acting on the unit cell
FiniteMPO identity_mpo(SiteListType const& SiteList);

// Returns the string MPO corresponding to the given local operator name
FiniteMPO string_mpo(SiteListType const& SiteList,
                     std::string const& OpName, QuantumNumbers::QuantumNumber const& Trans);

FiniteMPO string_mpo(SiteListType const& SiteList, std::string const& OpName);

// Given an expression, parse it as a SiteOperator on each site of SitesList,
// and form it into a scalar MPO.  The result is a string operator.
FiniteMPO
ParseStringOperator(SiteListType const& SiteList, std::string const& Expr, int Size);

// returns true if Op1 and Op2 are equal, to the specified tolerance
bool equal(FiniteMPO const& Op1, FiniteMPO const& Op2, double Tol = 1E-15);

// calculates the logarithm of the squared Frobenius norm of the operator
double
log_norm_frob_sq(FiniteMPO const& Op);

// returns the logarithm of the inner product <Op1|Op2> as
// <Op1|Op2> = Result.first * exp(Result.second)
// Result.first is a complex number on the unit circle.
std::pair<std::complex<double>, double>
log_inner_prod(FiniteMPO const& Op1, FiniteMPO const& Op2);

#include "finite_mpo.cc"

#endif
