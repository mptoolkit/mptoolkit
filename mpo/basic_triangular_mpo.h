// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/basic_triangular_mpo.h
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
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
//
// BasicTriangularMPO: a representation for lattice operators that are in upper triangular form.
//
// It is up to the user to ensure that the TriangularOperator stays in uper-triangular form.
// All functions defined in this header are OK though, the only way to generate a non-upper-triangular
// operator is to modify the components by hand (don't do that!).
//
// Some operations return a "1x1" MPO.  This is an MPO where the Basis1() and Basis2() both have
// dimension 1.  If the operator has a non-trivial unit cell then it may be that some internal
// dimensions are larger than 1.

#if !defined(MPTOOLKIT_MPO_BASIC_TRIANGULAR_MPO_H)
#define MPTOOLKIT_MPO_BASIC_TRIANGULAR_MPO_H

#include "generic_mpo.h"
#include "basic_finite_mpo.h"
#include <ostream>

class ProductMPO;

class BasicTriangularMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type;
      typedef data_type::iterator       iterator;
      typedef data_type::const_iterator const_iterator;
      typedef value_type::basis1_type   basis_type;

      BasicTriangularMPO() {}

      explicit BasicTriangularMPO(int Size) : Data_(Size) {}

      // construction as a single-site operator
      explicit BasicTriangularMPO(value_type const& x) : Data_(1, x) {}

      explicit BasicTriangularMPO(std::vector<value_type> const& x) : Data_(x.begin(), x.end()) {}

      int size() const { return Data_.size(); }

      bool empty() const { return Data_.empty(); }

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      // The Basis1() and Basis2() always coincide for a BasicTriangularMPO
      basis_type Basis() const { return Data_.front().Basis1(); }
      basis_type Basis1() const { return Data_.front().Basis1(); }
      basis_type Basis2() const { return Data_.back().Basis2(); }

      value_type const& operator[](int n) const { return Data_[n]; }
      value_type& operator[](int n) { return Data_[n]; }

      value_type const& front() const { return Data_.front(); }
      value_type& front() { return Data_.front(); }

      value_type const& back() const { return Data_.back(); }
      value_type& back() { return Data_.back(); }

      QuantumNumber TransformsAs() const { return this->Basis().front(); }

      // returns the component at entry (i,j).  Result is a 1x1 MPO (although the internal bond dimension might be bigger)
      BasicFiniteMPO operator()(int i, int j) const;

      // Returns the 1x1 MPO on the top left diagonal, the left 'string' term,
      // equivalent to operator()(0,0)
      BasicFiniteMPO left_string() const { return this->operator()(0,0); }

      // Returns the 1x1 MPO on the bottom right diagonal, the right 'string' term,
      // equivalent to operator()(Basis().size()-1, Basis().size())
      BasicFiniteMPO right_string() const { return this->operator()(this->Basis().size()-1, this->Basis().size()-1); }

      // Returns the 1x1 MPO at the top right element, which corresponds to the
      // value of the MPO with support within the unit cell
      BasicFiniteMPO as_finite() const { return this->operator()(0, this->Basis().size()-1); }

      operator GenericMPO const&() const { return Data_; }

      // extracts a diagonal operator.  This is either null, or a product of SimpleRedOperator's.
      std::vector<SimpleRedOperator> diagonal(int i) const;

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const { return Data_.LocalBasis1List(); }
      std::vector<BasisList> LocalBasis2List() const { return Data_.LocalBasis2List(); }

      data_type const& data() const { return Data_; }

      //data_type& data() { return Data_; }  // dangerous, probably shouldn't exist

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Data_.GetSymmetryList(); }

      void check_structure() const;
      void debug_check_structure() const;

   private:
      data_type Data_;

      friend PStream::opstream& operator<<(PStream::opstream& out, BasicTriangularMPO const& Op);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, BasicTriangularMPO& Op);
};

std::ostream&
operator<<(std::ostream& out, BasicTriangularMPO const& op);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(BasicTriangularMPO const& Op, std::ostream& out, double UnityEpsilon, int Verbose = 0);

inline
void print_structure(BasicTriangularMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

#if 0
// extracts a single column from a triangular operator.  Result is an Nx1 row-vector operator
GenericMPO extract_column(BasicTriangularMPO const& Op, int Col);

// extracts a single column from a triangular operator, excluding the diagonal.  Result is an Nx1 row-vector operator
GenericMPO extract_lower_column(BasicTriangularMPO const& Op, int Col);
#endif

// Helper function to make a list of identity operators over a unit cell
std::vector<SimpleOperator>
MakeIdentityUnitCell(std::vector<BasisList> const& Sites);

// replicate an operator on a unit cell this many times
BasicTriangularMPO repeat(BasicTriangularMPO const& x, int Count);

// Two triangular operators are *compatible* if they have the same operator on the
// top-left and bottom-right entries.
bool is_compatible(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

// Multiply by scalar
BasicTriangularMPO operator*(BasicTriangularMPO const& Op, double x);
BasicTriangularMPO operator*(BasicTriangularMPO const& Op, std::complex<double> x);
BasicTriangularMPO operator*(double x, BasicTriangularMPO const& Op);
BasicTriangularMPO operator*(std::complex<double> x, BasicTriangularMPO const& Op);

BasicTriangularMPO& operator*=(BasicTriangularMPO& Op, double x);
BasicTriangularMPO& operator*=(BasicTriangularMPO& Op, std::complex<double> x);

// Addition of triangular operators.  This is only possible if the operators
// are compatible.
BasicTriangularMPO& operator+=(BasicTriangularMPO& Op, BasicTriangularMPO const& x);
BasicTriangularMPO& operator-=(BasicTriangularMPO& Op, BasicTriangularMPO const& x);

BasicTriangularMPO operator+(BasicTriangularMPO const& x, BasicTriangularMPO const& y);
BasicTriangularMPO operator-(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

// unary negation
BasicTriangularMPO operator-(BasicTriangularMPO const& x);

// does a N-1 coarse-graining of the operator
BasicTriangularMPO coarse_grain(BasicTriangularMPO const& x, int N);

// Multiplication of triangular MPO's.  This doesn't depend on the
// compatibility of the operators.
BasicTriangularMPO prod(BasicTriangularMPO const& x, BasicTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);

// Disjoint prod.  Ignores the internal identity operators on the diagonal.
// Really only useful for commutators, since [A,B] = disjoint_prod(A,B) - disjoint_prod(B,A).
// If the MPO's are not first order then this will still leave additional identity components on the diagonal.
BasicTriangularMPO disjoint_prod(BasicTriangularMPO const& x, BasicTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);

// commutator x,y]
BasicTriangularMPO commutator(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
BasicTriangularMPO dot(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
BasicTriangularMPO inner(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

// cross product (if it exists)
BasicTriangularMPO cross(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
BasicTriangularMPO outer(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

BasicTriangularMPO pow(BasicTriangularMPO const& x, int n);

BasicTriangularMPO operator*(BasicTriangularMPO const& x, BasicTriangularMPO const& y);

BasicTriangularMPO& operator*=(BasicTriangularMPO& x, BasicTriangularMPO const& y);

// Get initial (1x1) E and F matrices.
StateComponent Initial_E(BasicTriangularMPO const& m);
StateComponent Initial_F(BasicTriangularMPO const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(BasicTriangularMPO const& m, VectorBasis const& B);
StateComponent Initial_F(BasicTriangularMPO const& m, VectorBasis const& B);

// Split an MPO into local operators.
// Result'[i] is a vector of 1x1 MPO's that have support starting from site i
// the length of the MPO is the number of sites where the operator is supported.
std::vector<std::vector<BasicFiniteMPO>>
SplitMPO(BasicTriangularMPO const& m);

// *Not yet implemented
// Analyze the dependencies of the MPO, for the E matrices.
// This returns a vector listing the largest row number (for E matrices) or
// smallest row number (for F matrices) that the given row/column depends on.
// For example, given the MPO
// (1 a)
// (0 1)
// the E-dependences are {-1,0},
// because row 0 of the E matrix doesn't have any dependencies, and
// calculating row 1 requires row 0.
// Similarly the F-dependences are {1,2}, since column 1 doesn't have
// any dependencies (column 2 is off the end), and column 0 depends on column 1.
// For an MPO with non-trivial diagonals, it is more interesting, eg
//
// (1 x y z a)
// (0 0 a b 0)
// (0 0 1 0 c)
// (0 0 0 1 d)
// (0 0 0 0 1)
// Here the E-dependencies are {-1, 0, 1, 1, 3}.  In particular, note that
// the second diagonal element doesn't require that the first diagonal element is
// already calculated, since it only uses E-elements [0,1].
std::vector<int> Dependencies_E(BasicTriangularMPO const& m);

BasicTriangularMPO gauge_flip(BasicTriangularMPO const& Op);

void optimize(BasicTriangularMPO& Op);

// optimize the representation using qr_decomposition
void qr_optimize(BasicTriangularMPO& Op);

// balances a triangular MPO - gives terms the same operator norm from
// the left and the right.
void balance(BasicTriangularMPO& Op);

// calculates the logarithm of the squared Frobenius norm of the operator
double
log_norm_frob_sq(BasicTriangularMPO const& Op);

// returns the logarithm of the inner product <Op1|Op2> as
// <Op1|Op2> = Result.first * exp(Result.second)
// Result.first is a complex number on the unit circle.
std::pair<std::complex<double>, double>
log_inner_prod(BasicTriangularMPO const& Op1, BasicTriangularMPO const& Op2);

// returns true if Op1 and Op2 are equal to the specified tolerance
bool
equal(BasicTriangularMPO const& Op1, BasicTriangularMPO const& Op2, double Tol);

inline
void
BasicTriangularMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

namespace LinearAlgebra
{

template <>
struct interface<BasicTriangularMPO>
{
   typedef void type;
};

template <>
struct Herm<BasicTriangularMPO>
{
   typedef BasicTriangularMPO const& argument_type;
   typedef HermitianProxy<BasicTriangularMPO> result_type;

   result_type operator()(argument_type x) const
   { return result_type(x); }
};

template <>
struct Conj<BasicTriangularMPO>
{
   typedef BasicTriangularMPO const& argument_type;
   typedef BasicTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      BasicTriangularMPO Result(x);
      for (BasicTriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = conj(*I);
      }
      return Result;
   }
};

template <>
struct Adjoint<BasicTriangularMPO>
{
   typedef BasicTriangularMPO const& argument_type;
   typedef BasicTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      BasicTriangularMPO Result(x);
      for (BasicTriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = adjoint(*I);
      }
      return Result;
   }
};

template <>
struct InvAdjoint<BasicTriangularMPO>
{
   typedef BasicTriangularMPO const& argument_type;
   typedef BasicTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      BasicTriangularMPO Result(x);
      for (BasicTriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = inv_adjoint(*I);
      }
      return Result;
   }
};

} // namespace LinearAlgebra

#endif
