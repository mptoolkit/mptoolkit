// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/infinite_mpo.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

// InfiniteMPO
//
// A variant of complex, ZeroMPO, BasicTriangularMPO and ProductMPO.
// The complex can also be interpreted as a constant multiplied by the
// identity ProductMPO.

#if !defined(MPTOOLKIT_MPO_INFINITE_MPO_H)
#define MPTOOLKIT_MPO_INFINITE_MPO_H

#include "zero_mpo.h"
#include "basic_triangular_mpo.h"
#include "product_mpo.h"
#include <boost/variant.hpp>

typedef boost::variant<std::complex<double>, ZeroMPO, BasicTriangularMPO, ProductMPO>
InfiniteMPOElement;

class InfiniteMPO
{
   public:
      typedef boost::variant<std::complex<double>, ZeroMPO, BasicTriangularMPO, ProductMPO> operator_type;

      InfiniteMPO() : Operator(ZeroMPO()) {}

      InfiniteMPO(operator_type const& Op) : Operator(Op) {}
      InfiniteMPO(InfiniteMPO const& Op) : Operator(Op.Operator), Description(Op.Description) {}
      InfiniteMPO(std::complex<double> const& x) : Operator(x) {}
      InfiniteMPO(double x) : Operator(std::complex<double>(x)) {}
      InfiniteMPO(BasicTriangularMPO const& Op) : Operator(Op) {}
      InfiniteMPO(ProductMPO const& Op) : Operator(Op) {}

      InfiniteMPO& operator=(InfiniteMPO const& Op) { Operator=Op.op(); Description = Op.Description; return *this; }
      InfiniteMPO& operator=(std::complex<double> const& x) { Operator=x; return *this; }
      InfiniteMPO& operator=(double x) { Operator=std::complex<double>(x); return *this; }
      InfiniteMPO& operator=(BasicTriangularMPO const& Op) { Operator=Op; return *this; }
      InfiniteMPO& operator=(ProductMPO const& Op) { Operator=Op; return *this; }

      // returns a friendly name for this object, either 'complex', 'BasicTriangularMPO',
      // 'ProductMPO'
      std::string name() const;

      // returns true if this operator is a BasicTriangularMPO
      bool is_triangular() const;

      BasicTriangularMPO const& as_basic_triangular_mpo() const;

      // returns true if this operator is a ProductMPO
      bool is_product() const;

      ProductMPO const& as_product_mpo() const;

      // returns true if this operator is a c-number
      bool is_complex() const;

      std::complex<double> as_complex() const;

      GenericMPO const& as_generic_mpo() const;

      operator_type& op() { return Operator; }
      operator_type const& op() const { return Operator; }


      void set_description(std::string const& s)
      { Description = s; }

      std::string const& description() const
      { return Description; }

#if 0
      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v) const
      {
         return Operator.apply_visitor(v);
      }

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v)
      {
         return Operator.apply_visitor(v);
      }
#endif

      static PStream::VersionTag VersionT;

   private:
      operator_type Operator;
      std::string Description;

      friend PStream::opstream& operator<<(PStream::opstream& out, InfiniteMPO const& Op);
      friend PStream::ipstream& operator>>(PStream::ipstream& in, InfiniteMPO& Op);
};

// operations.  The implementations can re-use the visitor patterns for the parser.
// None of these operations can mix BasicTriangularMPO and ProductMPO forms.

std::ostream& operator<<(std::ostream& out, InfiniteMPO const& Op);

InfiniteMPO& operator*=(InfiniteMPO& x, double a);
InfiniteMPO& operator*=(InfiniteMPO& x, std::complex<double> a);

InfiniteMPO& operator+=(InfiniteMPO& x, InfiniteMPO const& y);
InfiniteMPO& operator-=(InfiniteMPO& x, InfiniteMPO const& y);

InfiniteMPO operator+(InfiniteMPO const& x, InfiniteMPO const& y);
InfiniteMPO operator-(InfiniteMPO const& x, InfiniteMPO const& y);

InfiniteMPO operator-(InfiniteMPO const& x);

InfiniteMPO operator*(double a, InfiniteMPO const& x);
InfiniteMPO operator*(InfiniteMPO const& x, double a);
InfiniteMPO operator*(std::complex<double> a, InfiniteMPO const& x);
InfiniteMPO operator*(InfiniteMPO const& x, std::complex<double> a);

InfiniteMPO prod(InfiniteMPO const& x, InfiniteMPO const& y, QuantumNumbers::QuantumNumber const& q);
InfiniteMPO prod(InfiniteMPO const& x, InfiniteMPO const& y);
InfiniteMPO operator*(InfiniteMPO const& x, InfiniteMPO const& y);

InfiniteMPO commutator(InfiniteMPO const& x, InfiniteMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
InfiniteMPO dot(InfiniteMPO const& x, InfiniteMPO const& y);

// inner product - equivalent to dot(adjoint(x),y)
InfiniteMPO inner(InfiniteMPO const& x, InfiniteMPO const& y);

// cross product (if it exists)
InfiniteMPO cross(InfiniteMPO const& x, InfiniteMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
InfiniteMPO outer(InfiniteMPO const& x, InfiniteMPO const& y);

// power of an operator.  Defined for complex x and any n,
// or x BasicTriangularMPO and n >= 0
InfiniteMPO pow(InfiniteMPO const& x, int n);

// power of an operator.  Defined for complex x,y,
// or x BasicTriangularMPO and y integer >= 0
InfiniteMPO pow(InfiniteMPO const& x, InfiniteMPO const& y);

// Exponential - only defined for complex
InfiniteMPO exp(InfiniteMPO const& x);

// Conjugate
InfiniteMPO conj(InfiniteMPO const& x);

// Adjoint
InfiniteMPO adjoint(InfiniteMPO const& x);

// Inverse Adjoint
InfiniteMPO inv_adjoint(InfiniteMPO const& x);

#endif
