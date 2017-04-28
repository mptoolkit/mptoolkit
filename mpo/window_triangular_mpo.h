// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/window_triangular_mpo.h
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
// WindowTriangularMPO: a representation for lattice operators that are in upper triangular form.
//
// It is up to the user to ensure that the TriangularOperator stays in uper-triangular form.
// All functions defined in this header are OK though, the only way to generate a non-upper-triangular
// operator is to modify the components by hand (don't do that!).
//
// Some operations return a "1x1" MPO.  This is an MPO where the Basis1() and Basis2() both have
// dimension 1.  If the operator has a non-trivial unit cell then it may be that some internal
// dimensions are larger than 1.

#if !defined(MPTOOLKIT_MPO_WINDOW_TRIANGULAR_MPO_H)
#define MPTOOLKIT_MPO_WINDOW_TRIANGULAR_MPO_H

#include "basic_triangular_mpo.h"
#include "infinite_triangular_mpo.h"
#include "finite_triangular_mpo.h"
#include "finite_mpo.h"
#include <ostream>

class WindowTriangularMPO
{
   public:
      typedef OperatorComponent        value_type;

      WindowTriangularMPO() {}

      explicit WindowTriangularMPO(InfiniteTriangularMPO const& Other);

      int window_size() const { return Data_.size(); }

      BasicTriangularMPO const& left() { return Left; }
      BasicTriangularMPO const& window() { return Window; }
      BasicTriangularMPO const& right() { return Right; }

      bool empty() const { return Left.empty() && Window.empty() && Right.empty() }

      // The Basis1() and Basis2() are at the  left and  right edges of the window.
      // These coincide with the left and right semi-infinite MPO's respectively.
      basis_type Basis1() const { return Window.front().Basis1(); }
      basis_type Basis2() const { return Window.back().Basis2(); }

      // returns the site number of the first site of the window
      int window_offset_site() const { return WindowOffset; }

      // returns the unit cell number of the first unit cell of the window
      int window_offset_cell() const { return WindowOffset / Left.size(); }

      value_type const& operator[](int n) const { return Data_[n]; }
      value_type& operator[](int n) { return Data_[n]; }

      QuantumNumber TransformsAs() const { return this->Basis1().front(); }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const { return Data_.LocalBasis1List(); }
      std::vector<BasisList> LocalBasis2List() const { return Data_.LocalBasis2List(); }

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Data_.GetSymmetryList(); }

      void check_structure() const;
      void debug_check_structure() const;

      void set_components(BasicTriangularMPO const& L, BasicTriangularMPO const& W, BasicTriangularMPO const& R);

   private:
      // site index of the first site of the window
      int WindowOffset;
      BasicTriangularMPO Left;
      BasicTriangularMPO Window;
      BasicTriangularMPO Right;

      friend PStream::opstream& operator<<(PStream::opstream& out, WindowTriangularMPO const& Op);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, WindowTriangularMPO& Op);
};

std::ostream&
operator<<(std::ostream& out, WindowTriangularMPO const& op);

// prints the structure of the component, as an 'x' for a non-zero
// component or blank for a zero component
void print_structure(WindowTriangularMPO const& Op, std::ostream& out, double UnityEpsilon, int Verbose = 0);

inline
void print_structure(WindowTriangularMPO const& Op, std::ostream& out)
{
   print_structure(Op, out, DefaultClassifyUnityEpsilon);
}

// Two triangular operators are *compatible* if they have the same operator on the
// top-left and bottom-right entries.
bool is_compatible(WindowTriangularMPO const& x, WindowTriangularMPO const& y);

// Multiply by scalar
WindowTriangularMPO operator*(WindowTriangularMPO const& Op, double x);
WindowTriangularMPO operator*(WindowTriangularMPO const& Op, std::complex<double> x);
WindowTriangularMPO operator*(double x, WindowTriangularMPO const& Op);
WindowTriangularMPO operator*(std::complex<double> x, WindowTriangularMPO const& Op);

WindowTriangularMPO& operator*=(WindowTriangularMPO& Op, double x);
WindowTriangularMPO& operator*=(WindowTriangularMPO& Op, std::complex<double> x);

// Addition of triangular operators.  This is only possible if the operators
// are compatible.
WindowTriangularMPO& operator+=(WindowTriangularMPO& Op, WindowTriangularMPO const& x);
WindowTriangularMPO& operator-=(WindowTriangularMPO& Op, WindowTriangularMPO const& x);

WindowTriangularMPO& operator+=(WindowTriangularMPO& Op, FiniteTriangularMPO const& x);
WindowTriangularMPO& operator-=(WindowTriangularMPO& Op, FiniteTriangularMPO const& x);

WindowTriangularMPO& operator+=(WindowTriangularMPO& Op, BasicFiniteMPO const& x);
WindowTriangularMPO& operator-=(WindowTriangularMPO& Op, BasicFiniteMPO const& x);

WindowTriangularMPO operator+(WindowTriangularMPO const& x, WindowTriangularMPO const& y);
WindowTriangularMPO operator-(WindowTriangularMPO const& x, WindowTriangularMPO const& y);

WindowTriangularMPO operator+(WindowTriangularMPO const& x, BasicFiniteMPO const& y);
WindowTriangularMPO operator-(WindowTriangularMPO const& x, BasicFiniteMPO const& y);

WindowTriangularMPO operator+(BasicFiniteMPO const& x, WindowTriangularMPO const& y);
WindowTriangularMPO operator-(BasicFiniteMPO const& x, WindowTriangularMPO const& y);

// multiplication
WindowTriangularMPO operator*(WindowTriangularMPO const& x, WindowTriangularMPO const& y);
WindowTriangularMPO operator*(WindowTriangularMPO const& x, InfiniteTriangularMPO const& y);
WindowTriangularMPO operator*(InfiniteTriangularMPO const& x, WindowTriangularMPO const& y);
WindowTriangularMPO operator*(WindowTriangularMPO const& x, BasicFiniteMPO const& y);
WindowTriangularMPO operator*(BasicFiniteMPO const& x, WindowTriangularMPO const& y);
WindowTriangularMPO operator*(InfiniteTriangularMPO const& x, BasicFiniteMPO const& y);
WindowTriangularMPO operator*(BasicFiniteMPO const& x, InfiniteTriangularMPO const& y);

WindowTriangularMPO& operator*=(WindowTriangularMPO& x, WindowTriangularMPO const& y);
WindowTriangularMPO& operator*=(WindowTriangularMPO& x, InfiniteTriangularMPO const& y);
WindowTriangularMPO& operator*=(WindowTriangularMPO& x, BasicFiniteMPO const& y);

// TODO: we could define multiply by a WindowProductMPO
// TODO: do we also have a FiniteTriangularMPO ?

// unary negation
WindowTriangularMPO operator-(WindowTriangularMPO const& x);

// does a N-1 coarse-graining of the operator
WindowTriangularMPO coarse_grain(WindowTriangularMPO const& x, int N);

// Multiplication of triangular MPO's.
WindowTriangularMPO prod(WindowTriangularMPO const& x, WindowTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);
WindowTriangularMPO prod(WindowTriangularMPO const& x, InfiniteTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);
WindowTriangularMPO prod(InfiniteTriangularMPO const& x, WindowTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);
WindowTriangularMPO prod(WindowTriangularMPO const& x, BasicFiniteMPO const& y, QuantumNumbers::QuantumNumber const& q);
WindowTriangularMPO prod(BasicFiniteMPO const& x, WindowTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);
WindowTriangularMPO prod(InfiniteTriangularMPO const& x, BasicFiniteMPO const& y, QuantumNumbers::QuantumNumber const& q);
WindowTriangularMPO prod(BasicFiniteMPO const& x, InfiniteTriangularMPO const& y, QuantumNumbers::QuantumNumber const& q);

// dot product - takes into account the multiplicity to rescale the result
WindowTriangularMPO dot(WindowTriangularMPO const& x, WindowTriangularMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
WindowTriangularMPO inner(WindowTriangularMPO const& x, WindowTriangularMPO const& y);

// cross product (if it exists)
WindowTriangularMPO cross(WindowTriangularMPO const& x, WindowTriangularMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
WindowTriangularMPO outer(WindowTriangularMPO const& x, WindowTriangularMPO const& y);

WindowTriangularMPO pow(WindowTriangularMPO const& x, int n);

// Get initial (1x1) E and F matrices.
StateComponent Initial_E(WindowTriangularMPO const& m);
StateComponent Initial_F(WindowTriangularMPO const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(WindowTriangularMPO const& m, VectorBasis const& B);
StateComponent Initial_F(WindowTriangularMPO const& m, VectorBasis const& B);

void optimize(WindowTriangularMPO& Op);

// calculates the logarithm of the squared Frobenius norm of the operator
double
log_norm_frob_sq(WindowTriangularMPO const& Op);

// returns the logarithm of the inner product <Op1|Op2> as
// <Op1|Op2> = Result.first * exp(Result.second)
// Result.first is a complex number on the unit circle.
std::pair<std::complex<double>, double>
log_inner_prod(WindowTriangularMPO const& Op1, WindowTriangularMPO const& Op2);

// returns true if Op1 and Op2 are equal to the specified tolerance
bool
equal(WindowTriangularMPO const& Op1, WindowTriangularMPO const& Op2, double Tol);

inline
void
WindowTriangularMPO::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

namespace LinearAlgebra
{

template <>
struct interface<WindowTriangularMPO>
{
   typedef void type;
};

template <>
struct Herm<WindowTriangularMPO>
{
   typedef WindowTriangularMPO const& argument_type;
   typedef HermitianProxy<WindowTriangularMPO> result_type;

   result_type operator()(argument_type x) const
   { return result_type(x); }
};

template <>
struct Conj<WindowTriangularMPO>
{
   typedef WindowTriangularMPO const& argument_type;
   typedef WindowTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      WindowTriangularMPO Result(x);
      Result.set_components(conj(x.left()), conj(x.window()), conj(x.right()));
      return Result;
   }
};

template <>
struct Adjoint<WindowTriangularMPO>
{
   typedef WindowTriangularMPO const& argument_type;
   typedef WindowTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      WindowTriangularMPO Result(x);
      Result.set_components(adjoint(x.left()), adjoint(x.window()), adjoint(x.right()));
      return Result;
   }
};

template <>
struct InvAdjoint<WindowTriangularMPO>
{
   typedef WindowTriangularMPO const& argument_type;
   typedef WindowTriangularMPO result_type;

   result_type operator()(argument_type x) const
   {
      WindowTriangularMPO Result(x);
      Result.set_components(inv_adjoint(x.left()), inv_adjoint(x.window()), inv_adjoint(x.right()));
      return Result;
   }
};

} // namespace LinearAlgebra

#endif
