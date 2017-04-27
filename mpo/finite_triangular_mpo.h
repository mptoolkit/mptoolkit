// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/triangular_mpo.h
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

#if !defined(MPTOOLKIT_MPO_FINITE_TRIANGULAR_MPO_H)
#define MPTOOLKIT_MPO_FINITE_TRIANGULAR_MPO_H

#include "basic_triangular_mpo.h"

class FiniteTriangularMPO
{
   private:
      typedef BasicTriangularMPO data_type;

   public:
      typedef OperatorComponent         value_type;
      typedef data_type::iterator       iterator;
      typedef data_type::const_iterator const_iterator;
      typedef value_type::basis1_type   basis_type;

      FiniteTriangularMPO() {}

      explicit FiniteTriangularMPO(int Size) : Data_(Size) {}

      // construction as a single-site operator
      explicit FiniteTriangularMPO(value_type const& x) : Data_(1, x) {}

      explicit FiniteTriangularMPO(std::vector<value_type> const& x) : Data_(x.begin(), x.end()) {}

      int size() const { return Data_.size(); }

      bool empty() const { return Data_.empty(); }

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      basis_type Basis1() const { return Data_.front().Basis1(); }
      basis_type Basis2() const { return Data_.back().Basis2(); }

      value_type const& operator[](int n) const { return Data_[n]; }
      value_type& operator[](int n) { return Data_[n]; }

      value_type const& front() const { return Data_.front(); }
      value_type& front() { return Data_.front(); }

      value_type const& back() const { return Data_.back(); }
      value_type& back() { return Data_.back(); }

      QuantumNumber TransformsAs() const { return this->Basis1().front(); }

      operator BasicTriangularMPO const&() const { return Data_; }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const { return Data_.LocalBasis1List(); }
      std::vector<BasisList> LocalBasis2List() const { return Data_.LocalBasis2List(); }

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Data_.GetSymmetryList(); }

      void check_structure() const;
      void debug_check_structure() const;

   private:
      data_type Data_;

      friend PStream::opstream& operator<<(PStream::opstream& out, FiniteTriangularMPO const& Op);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, FiniteTriangularMPO& Op);
};

std::ostream&
operator<<(std::ostream& out, FiniteTriangularMPO const& op);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(FiniteTriangularMPO const& Op, std::ostream& out, double UnityEpsilon, int Verbose = 0);

inline
void print_structure(FiniteTriangularMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

// replicate an operator on a unit cell this many times
FiniteTriangularMPO repeat(FiniteTriangularMPO const& x, int Count);

// Two triangular operators are *compatible* if they have the same operator on the
// top-left and bottom-right entries.
bool is_compatible(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

// Multiply by scalar
FiniteTriangularMPO operator*(FiniteTriangularMPO const& Op, double x);
FiniteTriangularMPO operator*(FiniteTriangularMPO const& Op, std::complex<double> x);
FiniteTriangularMPO operator*(double x, FiniteTriangularMPO const& Op);
FiniteTriangularMPO operator*(std::complex<double> x, FiniteTriangularMPO const& Op);

FiniteTriangularMPO& operator*=(FiniteTriangularMPO& Op, double x);
FiniteTriangularMPO& operator*=(FiniteTriangularMPO& Op, std::complex<double> x);

// Addition of triangular operators.  This is only possible if the operators
// are compatible.
FiniteTriangularMPO& operator+=(FiniteTriangularMPO& Op, FiniteTriangularMPO const& x);
FiniteTriangularMPO& operator-=(FiniteTriangularMPO& Op, FiniteTriangularMPO const& x);

FiniteTriangularMPO operator+(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);
FiniteTriangularMPO operator-(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

// unary negation
FiniteTriangularMPO operator-(FiniteTriangularMPO const& x);

// does a N-1 coarse-graining of the operator
FiniteTriangularMPO coarse_grain(FiniteTriangularMPO const& x, int N);

// Multiplication of triangular MPO's.  This doesn't depend on the
// compatibility of the operators.
FiniteTriangularMPO prod(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);

// dot product - takes into account the multiplicity to rescale the result
FiniteTriangularMPO dot(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
FiniteTriangularMPO inner(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

// cross product (if it exists)
FiniteTriangularMPO cross(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
FiniteTriangularMPO outer(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

FiniteTriangularMPO pow(FiniteTriangularMPO const& x, int n);

FiniteTriangularMPO operator*(FiniteTriangularMPO const& x, FiniteTriangularMPO const& y);

FiniteTriangularMPO& operator*=(FiniteTriangularMPO& x, FiniteTriangularMPO const& y);

// Get initial (1x1) E and F matrices.
StateComponent Initial_E(FiniteTriangularMPO const& m);
StateComponent Initial_F(FiniteTriangularMPO const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(FiniteTriangularMPO const& m, VectorBasis const& B);
StateComponent Initial_F(FiniteTriangularMPO const& m, VectorBasis const& B);

void optimize(FiniteTriangularMPO& Op);

// calculates the logarithm of the squared Frobenius norm of the operator
double
log_norm_frob_sq(FiniteTriangularMPO const& Op);

// returns the logarithm of the inner product <Op1|Op2> as
// <Op1|Op2> = Result.first * exp(Result.second)
// Result.first is a complex number on the unit circle.
std::pair<std::complex<double>, double>
log_inner_prod(FiniteTriangularMPO const& Op1, FiniteTriangularMPO const& Op2);

// returns true if Op1 and Op2 are equal to the specified tolerance
bool
equal(FiniteTriangularMPO const& Op1, FiniteTriangularMPO const& Op2, double Tol);

inline
void
FiniteTriangularMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

namespace LinearAlgebra
{

template <>
struct interface<FiniteTriangularMPO>
{
   typedef void type;
};

template <>
struct Herm<FiniteTriangularMPO>
{
   typedef FiniteTriangularMPO const& argument_type;
   typedef HermitianProxy<FiniteTriangularMPO> result_type;

   result_type operator()(argument_type x) const
   { return result_type(x); }
};

template <>
struct Conj<FiniteTriangularMPO>
{
   typedef FiniteTriangularMPO const& argument_type;
   typedef FiniteTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      FiniteTriangularMPO Result(x);
      for (FiniteTriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = conj(*I);
      }
      return Result;
   }
};

template <>
struct Adjoint<FiniteTriangularMPO>
{
   typedef FiniteTriangularMPO const& argument_type;
   typedef FiniteTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      FiniteTriangularMPO Result(x);
      for (FiniteTriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = adjoint(*I);
      }
      return Result;
   }
};

template <>
struct InvAdjoint<FiniteTriangularMPO>
{
   typedef FiniteTriangularMPO const& argument_type;
   typedef FiniteTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      FiniteTriangularMPO Result(x);
      for (FiniteTriangularMPO::iterator I = Result.begin(); I != Result.end(); ++I)
      {
         *I = inv_adjoint(*I);
      }
      return Result;
   }
};

} // namespace LinearAlgebra

#endif
