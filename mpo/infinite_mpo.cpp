// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/infinite_mpo.cpp
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "infinite_mpo.h"
#include "parser/visitor_actions.h"
#include "mpo/infinite_mpo_actions.h"
#include "pstream/variant.h"

// Warning: Due to the implicit conversion of BasicTriangularMPO and ProductMPO to
// InfiniteMPOElement and InfiniteMPO, we can get some weird errors if
// the visitor functions don't find a match.

// Streaming versions:
//
// No explicit versioning information, check the LatticeVersion for <= 4
// std::variant<std::complex<double>, BasicTriangularMPO, ProductMPO> Operator
// std::string                                                    Description
//
// Version 2:
// Explicit version support.
// std::variant<std::complex<double>, ZeroMPO, BasicTriangularMPO, ProductMPO> Operator
// std::string                                                             Description

PStream::VersionTag
InfiniteMPO::VersionT(2);

using namespace Parser;

PStream::opstream&
operator<<(PStream::opstream& out, InfiniteMPO const& Op)
{
   out << InfiniteMPO::VersionT.default_version();
   out << Op.Operator;
   out << Op.Description;
   return out;
}

extern PStream::VersionTag LatticeVersion;

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteMPO& Op)
{
   // Prior to version 5, we didn't have separate versioning support for InfiniteMPO.
   int Version = in.version_of(LatticeVersion) >= 5 ? in.read<int>() : 1;

   if (Version == 1)
   {
      std::variant<std::complex<double>, BasicTriangularMPO, ProductMPO> x;
      in >> x;
      std::visit([&Op](auto const& Value) { Op.Operator = Value; }, x);
      in >> Op.Description;
      return in;
   }
   // else

   in >> Op.Operator;
   in >> Op.Description;
   return in;
}

std::string
InfiniteMPO::name() const
{
   if (this->is_complex())
      return "complex";
   if (this->is_triangular())
      return "BasicTriangularMPO";
   if (this->is_product())
      return "ProductMPO";
   return "unknown";
}

bool
InfiniteMPO::is_triangular() const
{
   return bool(std::get_if<BasicTriangularMPO>(&this->op()));
}

BasicTriangularMPO const&
InfiniteMPO::as_basic_triangular_mpo() const
{
   return std::get<BasicTriangularMPO>(this->op());
}

bool
InfiniteMPO::is_product() const
{
   return bool(std::get_if<ProductMPO>(&this->op()));
}

ProductMPO const&
InfiniteMPO::as_product_mpo() const
{
   return std::get<ProductMPO>(this->op());
}

bool
InfiniteMPO::is_complex() const
{
   return bool(std::get_if<std::complex<double> >(&this->op()));
}

std::complex<double>
InfiniteMPO::as_complex() const
{
   return std::get<std::complex<double> >(this->op());
}

GenericMPO const&
InfiniteMPO::as_generic_mpo() const
{
   if (this->is_triangular())
      return this->as_basic_triangular_mpo();
   else if (this->is_product())
      return this->as_product_mpo();
   else
      throw std::runtime_error("fatal: InfiniteMPO must be a triangular or product MPO to convert to GenericMPO.");
}

std::ostream& operator<<(std::ostream& out, InfiniteMPO const& Op)
{
   out << Op.name();
   return out;
}

InfiniteMPO& operator*=(InfiniteMPO& x, double a)
{
   x = x*a;
   return x;
}

InfiniteMPO& operator*=(InfiniteMPO& x, std::complex<double> a)
{
   x = x*a;
   return x;
}

InfiniteMPO& operator+=(InfiniteMPO& x, InfiniteMPO const& y)
{
   x = x+y;
   return x;
}

InfiniteMPO& operator-=(InfiniteMPO& x, InfiniteMPO const& y)
{
   x = x-y;
   return x;
}

InfiniteMPO operator+(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_addition<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO operator-(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_subtraction<InfiniteMPOElement>(), x.op(), y.op());
}


InfiniteMPO operator-(InfiniteMPO const& x)
{
   return std::visit(negate_element<InfiniteMPOElement>(), x.op());
}

InfiniteMPO operator*(double a, InfiniteMPO const& x)
{
   return InfiniteMPO(a)*x;
}

InfiniteMPO operator*(InfiniteMPO const& x, double a)
{
   return x*InfiniteMPO(a);
}

InfiniteMPO operator*(std::complex<double> a, InfiniteMPO const& x)
{
   return InfiniteMPO(a)*x;
}

InfiniteMPO operator*(InfiniteMPO const& x, std::complex<double> a)
{
   return x*InfiniteMPO(a);
}

InfiniteMPO prod(InfiniteMPO const& x, InfiniteMPO const& y,
                 QuantumNumbers::QuantumNumber const& q)
{
   return std::visit(ternary_product_q<InfiniteMPOElement>(q), x.op(), y.op());
}

InfiniteMPO prod(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_multiplication<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO operator*(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_multiplication<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO commutator(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_commutator<InfiniteMPOElement>(), x.op(), y.op());
}

// dot product - takes into account the multiplicity to rescale the result
InfiniteMPO dot(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_dot_product<InfiniteMPOElement>(), x.op(), y.op());
}

// inner product - equivalent to dot(adjoint(x),y)
InfiniteMPO inner(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_inner_product<InfiniteMPOElement>(), x.op(), y.op());
}

// cross product (if it exists)
InfiniteMPO cross(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_cross_product<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO outer(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_outer_product<InfiniteMPOElement>(), x.op(), y.op());
}

// power of an operator.  Requires n > 1.
InfiniteMPO pow(InfiniteMPO const& x, int n)
{
   return std::visit(unary_power<InfiniteMPOElement>(n), x.op());
}

InfiniteMPO pow(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return std::visit(binary_power<InfiniteMPOElement>(), x.op(), y.op());
}

// Exponential - only defined for complex
InfiniteMPO exp(InfiniteMPO const& x)
{
   return std::visit(ElementExp<InfiniteMPOElement>(), x.op());
}

// Conjugate
InfiniteMPO conj(InfiniteMPO const& x)
{
   return std::visit(ElementConj<InfiniteMPOElement>(), x.op());
}

// Adjoint
InfiniteMPO adjoint(InfiniteMPO const& x)
{
   return std::visit(ElementAdjoint<InfiniteMPOElement>(), x.op());
}

// Inverse Adjoint
InfiniteMPO inv_adjoint(InfiniteMPO const& x)
{
   return std::visit(ElementInvAdjoint<InfiniteMPOElement>(), x.op());
}
