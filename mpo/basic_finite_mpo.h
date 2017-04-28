// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/basic_finite_mpo.h
//
// Copyright (C) 2013-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
// generic operators.

#if !defined(MPTOOLKIT_MPO_FINITE_MPO_H)
#define MPTOOLKIT_MPO_FINITE_MPO_H

#include "generic_mpo.h"
#include "lattice/latticesite.h"

class BasicFiniteMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type;
      typedef data_type::const_iterator const_iterator;
      typedef data_type::iterator       iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      BasicFiniteMPO() {}

      BasicFiniteMPO(BasicFiniteMPO const& Other) : Data(Other.Data) {}

      // Removed this constructor because it doesn't make much sense to define a BasicFiniteMPO
      // without specifying the LatticeCommute
      //      explicit BasicFiniteMPO(int Size) : Data(Size) {}

      explicit BasicFiniteMPO(int Size) : Data(Size) {}

      // Construction from a generic MPO.  The generic MPO must already be in finite form.
      explicit BasicFiniteMPO(GenericMPO const& Other);

      BasicFiniteMPO& operator=(BasicFiniteMPO const& Other) { Data = Other.Data; return *this; }

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
      static BasicFiniteMPO make_identity(std::vector<BasisList> const& Basis);
      static BasicFiniteMPO make_identity(std::vector<BasisList> const& Basis,
                                     QuantumNumber const& q);

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, BasicFiniteMPO const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, BasicFiniteMPO& op);

// Returns the MPO that is Op1 \otimes Op2.
// PRECONDITION: Op1.Basis2() == Op2.Basis1()
BasicFiniteMPO join(BasicFiniteMPO const& Op1, BasicFiniteMPO const& Op2);

// Repeats Op Count number of times, Op^{\oprod Count}.
// PRECONDITION: Op.Basis2() == Op.Basis1()
BasicFiniteMPO repeat(BasicFiniteMPO const& Op, int Count);

BasicFiniteMPO& operator*=(BasicFiniteMPO& x, double a);
BasicFiniteMPO& operator*=(BasicFiniteMPO& x, std::complex<double> a);

BasicFiniteMPO& operator+=(BasicFiniteMPO& x, BasicFiniteMPO const& y);
BasicFiniteMPO& operator-=(BasicFiniteMPO& x, BasicFiniteMPO const& y);

BasicFiniteMPO operator+(BasicFiniteMPO const& x, BasicFiniteMPO const& y);
BasicFiniteMPO operator-(BasicFiniteMPO const& x, BasicFiniteMPO const& y);

BasicFiniteMPO operator-(BasicFiniteMPO const& x);

BasicFiniteMPO operator*(double a, BasicFiniteMPO const& x);
BasicFiniteMPO operator*(BasicFiniteMPO const& x, double a);
BasicFiniteMPO operator*(std::complex<double> a, BasicFiniteMPO const& x);
BasicFiniteMPO operator*(BasicFiniteMPO const& x, std::complex<double> a);

BasicFiniteMPO prod(BasicFiniteMPO const& x, BasicFiniteMPO const& y, QuantumNumbers::QuantumNumber const& q);
BasicFiniteMPO prod(BasicFiniteMPO const& x, BasicFiniteMPO const& y);
BasicFiniteMPO operator*(BasicFiniteMPO const& x, BasicFiniteMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
BasicFiniteMPO dot(BasicFiniteMPO const& x, BasicFiniteMPO const& y);

// cross product (if it exists)
BasicFiniteMPO cross(BasicFiniteMPO const& x, BasicFiniteMPO const& y);

// Helper function for the coupling coefficient for the outer() function
double outer_coefficient(int degree_x, int degree_y, int degree_q);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
BasicFiniteMPO outer(BasicFiniteMPO const& x, BasicFiniteMPO const& y);

// project a (reducible) operator onto an irreducible component
BasicFiniteMPO project(BasicFiniteMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.
BasicFiniteMPO pow(BasicFiniteMPO const& x, int n);

// Exponential operator.
BasicFiniteMPO exp(BasicFiniteMPO const& x);

// Conjugate
BasicFiniteMPO conj(BasicFiniteMPO const& x);

// Adjoint
BasicFiniteMPO adjoint(BasicFiniteMPO const& x);

// Inverse Adjoint
BasicFiniteMPO inv_adjoint(BasicFiniteMPO const& x);

// optimize the representation
void optimize(BasicFiniteMPO& Op);

// optimize the representation using QR decomposition
void qr_optimize(BasicFiniteMPO& Op);

// completely coarse-grain the MPO into a simple operator.
// The dimensions of this operator are exponentially big in the number of sites
// in x, so be careful!
// For non-abelian symmetries, this coarse-grain occurs from left to right.
SimpleRedOperator coarse_grain(BasicFiniteMPO const& x);

// Does a N-1 coarse graining of an operator.  The length must be a multiple of N
BasicFiniteMPO coarse_grain(BasicFiniteMPO const& Op, int N);

// The opposite of coarse_grain - decompose an operator acting on the entire Hilbert space
// into a BasicFiniteMPO
BasicFiniteMPO fine_grain(SimpleOperator const& x,
                     std::vector<BasisList> const& LocalBasis1,
                     std::vector<BasisList> const& LocalBasis2);

// Make an identity operator that acts on the same local Hilbert space as x
BasicFiniteMPO
MakeIdentityFrom(BasicFiniteMPO const& x);

// Make an identity operator that acts on the same local Hilbert space as x,
// with the given quantum number in the auxiliary basis
BasicFiniteMPO
MakeIdentityFrom(BasicFiniteMPO const& x, QuantumNumber const& q);

// output to a stream
std::ostream& operator<<(std::ostream& out, BasicFiniteMPO const& x);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(BasicFiniteMPO const& Op, std::ostream& out, double UnityEpsilon);

inline
void print_structure(BasicFiniteMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

// returns the BasicFiniteMPO for the identity operator acting on the unit cell
// with quantum number q in the auxiliary basis
BasicFiniteMPO identity_mpo(SiteListType const& SiteList, QuantumNumbers::QuantumNumber const& q);

// returns the BasicFiniteMPO for the identity operator acting on the unit cell
BasicFiniteMPO identity_mpo(SiteListType const& SiteList);

// Returns the string MPO corresponding to the given local operator name
BasicFiniteMPO string_mpo(SiteListType const& SiteList,
                     std::string const& OpName, QuantumNumbers::QuantumNumber const& Trans);

BasicFiniteMPO string_mpo(SiteListType const& SiteList, std::string const& OpName);

// Given an expression, parse it as a SiteOperator on each site of SitesList,
// and form it into a scalar MPO.  The result is a string operator.
BasicFiniteMPO
ParseStringOperator(SiteListType const& SiteList, std::string const& Expr, int Size);

// returns true if Op1 and Op2 are equal, to the specified tolerance
bool equal(BasicFiniteMPO const& Op1, BasicFiniteMPO const& Op2, double Tol = 1E-15);

// calculates the logarithm of the squared Frobenius norm of the operator
double
log_norm_frob_sq(BasicFiniteMPO const& Op);

// returns the logarithm of the inner product <Op1|Op2> as
// <Op1|Op2> = Result.first * exp(Result.second)
// Result.first is a complex number on the unit circle.
std::pair<std::complex<double>, double>
log_inner_prod(BasicFiniteMPO const& Op1, BasicFiniteMPO const& Op2);

#include "basic_finite_mpo.cc"

#endif
