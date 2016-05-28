// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/triangular_mpo.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER
//
// TriangularMPO: a representation for lattice operators that are in upper triangular form.
//
// It is up to the user to ensure that the TriangularOperator stays in uper-triangular form.
// All functions defined in this header are OK though, the only way to generate a non-upper-triangular
// operator is to modify the components by hand (don't do that!).
//
// Some operations return a "1x1" MPO.  This is an MPO where the Basis1() and Basis2() both have
// dimension 1.  If the operator has a non-trivial unit cell then it may be that some internal
// dimensions are larger than 1.

#if !defined(TRIANGULAR_MPO_H_SDJFHU47894789H9O43)
#define TRIANGULAR_MPO_H_SDJFHU47894789H9O43

#include "generic_mpo.h"
#include "finite_mpo.h"
#include <ostream>

class TriangularMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type;
      typedef data_type::iterator       iterator;
      typedef data_type::const_iterator const_iterator;
      typedef value_type::basis1_type   basis_type;

      TriangularMPO() {}

      explicit TriangularMPO(int Size) : Data_(Size) {}

      // construction as a single-site operator
      explicit TriangularMPO(value_type const& x) : Data_(1, x) {}

      explicit TriangularMPO(std::vector<value_type> const& x) : Data_(x.begin(), x.end()) {}

      int size() const { return Data_.size(); }

      bool empty() const { return Data_.empty(); }

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }
   
      // The Basis1() and Basis2() always coincide for a TriangularMPO
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
      FiniteMPO operator()(int i, int j) const;

      // Returns the 1x1 MPO on the top left diagonal, the left 'string' term,
      // equivalent to operator()(0,0)
      FiniteMPO left_string() const { return this->operator()(0,0); }

      // Returns the 1x1 MPO on the bottom right diagonal, the right 'string' term,
      // equivalent to operator()(Basis().size()-1, Basis().size())
      FiniteMPO right_string() const { return this->operator()(this->Basis().size()-1, this->Basis().size()-1); }

      // Returns the 1x1 MPO at the top right element, which corresponds to the
      // value of the MPO with support within the unit cell
      FiniteMPO as_finite() const { return this->operator()(0, this->Basis().size()-1); }

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

      friend PStream::opstream& operator<<(PStream::opstream& out, TriangularMPO const& Op);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, TriangularMPO& Op);
};

std::ostream&
operator<<(std::ostream& out, TriangularMPO const& op);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(TriangularMPO const& Op, std::ostream& out, double UnityEpsilon);

inline
void print_structure(TriangularMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

#if 0
// extracts a single column from a triangular operator.  Result is an Nx1 row-vector operator
GenericMPO extract_column(TriangularMPO const& Op, int Col);

// extracts a single column from a triangular operator, excluding the diagonal.  Result is an Nx1 row-vector operator
GenericMPO extract_lower_column(TriangularMPO const& Op, int Col);
#endif

TriangularMPO TriangularOneSite(SimpleOperator const& x);

// A one-site operator with the given momentum, in angular units
TriangularMPO TriangularOneSite(SimpleOperator const& x, double Momentum);

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y,
                                     QuantumNumbers::QuantumNumber const& Trans);

TriangularMPO TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y);

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, 
                                                std::complex<double> Factor, QuantumNumber const& Trans);

TriangularMPO TriangularTwoSiteExponential(SimpleOperator const& x, SimpleOperator const& y, 
                                                std::complex<double> Factor);

TriangularMPO TriangularThreeSite(SimpleOperator const& x, 
                                       SimpleOperator const& y, 
                                       SimpleOperator const& z);

// a two point interaction between sites n1 and n2 of a lattice.
// This probably should use SiteOperator's rather than SimpleOperator's, so that
// we can handle fermions, especially for the case n2<n1.
// this is assumed to be a bosonic operator; if n2 < n1 we swap the sites
// (which is incorrect for fermions)
TriangularMPO TwoPointOperator(std::vector<BasisList> const& Sites, 
                                    int n1, SimpleOperator const& x1,
                                    int n2, SimpleOperator const& x2);

// This is a restricted implementation; we must have
// PRECONDITION: 0 <= n1 < Sites.size()
// PRECONDITION: Sites.size() <= n2 < 2*Sites.size()
// That is, the operator1 must be in the first unit cell, the operator2 must be
// in the second unit cell.
TriangularMPO TwoPointExponentialOperator(std::vector<BasisList> const& Sites, 
                                          int n1, SimpleOperator const& x1,
                                          int n2, SimpleOperator const& x2,
                                          std::complex<double> Factor);


// a two-point string operator where the String term is inserted
// at sites n1+1, n1+2, ..., n2-1.
// Because this function can be used to implement fermionic operators,
// we demand normal ordering of sites; it is an error to call this function
// with n2<n1.
TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
                                     int n1, SimpleOperator const& x1,
                                     SimpleOperator const& String,
                                     int n2, SimpleOperator const& x2);

TriangularMPO TwoPointStringOperator(std::vector<BasisList> const& Sites, 
                                     int n1, SimpleOperator const& x1,
                                     std::vector<SimpleOperator> const& String,
                                     int n2, SimpleOperator const& x2);

// A one-site operator on a lattice with a given momentum, in angular units per unit cell
TriangularMPO OnePointOperator(std::vector<BasisList> const& Sites, 
                                    int n, SimpleOperator const& x, double Momentum = 0);

// A one-site operator on a lattice with a given momentum, in angular units per unit cell
TriangularMPO OnePointStringOperator(std::vector<BasisList> const& Sites, 
					  std::vector<SimpleOperator> const& String,
					  int n, SimpleOperator const& x, double Momentum = 0);

// A one-site operator on a lattice with a given momentum, in angular units per unit cell
// Here we exclude the string term on the unit cell where the operator is defined
TriangularMPO OnePointStringOperatorExclude(std::vector<BasisList> const& Sites, 
					    std::vector<SimpleOperator> const& String,
					    int n, SimpleOperator const& x, double Momentum = 0);

TriangularMPO OneCellStringOperator(std::vector<BasisList> const& Sites, 
				    std::vector<SimpleOperator> const& String,
				    std::vector<SimpleOperator> const& CellOp,
				    double Momentum = 0);

// Helper function to make a list of identity operators over a unit cell
std::vector<SimpleOperator>
MakeIdentityUnitCell(std::vector<BasisList> const& Sites);

// replicate an operator on a unit cell this many times
TriangularMPO repeat(TriangularMPO const& x, int Count);

// Two triangular operators are *compatible* if they have the same operator on the
// top-left and bottom-right entries.
bool is_compatible(TriangularMPO const& x, TriangularMPO const& y);

// Multiply by scalar
TriangularMPO operator*(TriangularMPO const& Op, double x);
TriangularMPO operator*(TriangularMPO const& Op, std::complex<double> x);
TriangularMPO operator*(double x, TriangularMPO const& Op);
TriangularMPO operator*(std::complex<double> x, TriangularMPO const& Op);

TriangularMPO& operator*=(TriangularMPO& Op, double x);
TriangularMPO& operator*=(TriangularMPO& Op, std::complex<double> x);

// Addition of triangular operators.  This is only possible if the operators
// are compatible.
TriangularMPO& operator+=(TriangularMPO& Op, TriangularMPO const& x);
TriangularMPO& operator-=(TriangularMPO& Op, TriangularMPO const& x);

TriangularMPO operator+(TriangularMPO const& x, TriangularMPO const& y);
TriangularMPO operator-(TriangularMPO const& x, TriangularMPO const& y);

// unary negation
TriangularMPO operator-(TriangularMPO const& x);

// does a N-1 coarse-graining of the operator
TriangularMPO coarse_grain(TriangularMPO const& x, int N);

// Multiplication of triangular MPO's.  This doesn't depend on the
// compatibility of the operators.
TriangularMPO prod(TriangularMPO const& x, TriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);

// dot product - takes into account the multiplicity to rescale the result
TriangularMPO dot(TriangularMPO const& x, TriangularMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
TriangularMPO inner(TriangularMPO const& x, TriangularMPO const& y);

// cross product (if it exists)
TriangularMPO cross(TriangularMPO const& x, TriangularMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
TriangularMPO outer(TriangularMPO const& x, TriangularMPO const& y);

TriangularMPO pow(TriangularMPO const& x, int n);

TriangularMPO operator*(TriangularMPO const& x, TriangularMPO const& y);

TriangularMPO& operator*=(TriangularMPO& x, TriangularMPO const& y);

// Get initial (1x1) E and F matrices.
StateComponent Initial_E(TriangularMPO const& m);
StateComponent Initial_F(TriangularMPO const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(TriangularMPO const& m, VectorBasis const& B);
StateComponent Initial_F(TriangularMPO const& m, VectorBasis const& B);

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
std::vector<int> Dependencies_E(TriangularMPO const& m);

void optimize(TriangularMPO& Op);

// balances a triangular MPO - gives terms the same operator norm from
// the left and the right.
void balance(TriangularMPO& Op);

inline
void
TriangularMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

namespace LinearAlgebra
{

template <>
struct interface<TriangularMPO>
{
   typedef void type;
};

template <>
struct Herm<TriangularMPO>
{
   typedef TriangularMPO const& argument_type;
   typedef HermitianProxy<TriangularMPO> result_type;

   result_type operator()(argument_type x) const 
   { return result_type(x); }
};

template <>
struct Conj<TriangularMPO>
{
   typedef TriangularMPO const& argument_type;
   typedef TriangularMPO result_type;

   result_type operator()(argument_type x) const 
   { 
      TriangularMPO Result(x);
      for (TriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
	 *I = conj(*I);
      }
      return Result;
   }
};

template <>
struct Adjoint<TriangularMPO>
{
   typedef TriangularMPO const& argument_type;
   typedef TriangularMPO result_type;

   result_type operator()(argument_type x) const 
   { 
      TriangularMPO Result(x);
      for (TriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
	 *I = adjoint(*I);
      }
      return Result;
   }
};

template <>
struct InvAdjoint<TriangularMPO>
{
   typedef TriangularMPO const& argument_type;
   typedef TriangularMPO result_type;

   result_type operator()(argument_type x) const 
   {
      TriangularMPO Result(x);
      for (TriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
	 *I = inv_adjoint(*I);
      }
      return Result;
   }
};

} // namespace LinearAlgebra

#endif
