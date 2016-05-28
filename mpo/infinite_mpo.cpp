// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/infinite_mpo.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "infinite_mpo.h"
#include "parser/visitor_actions.h"
#include "mpo/infinite_mpo_actions.h"
#include "pstream/variant.h"

// Warning: Due to the implicit conversion of TriangularMPO and ProductMPO to 
// InfiniteMPOElement and InfiniteMPO, we can get some weird errors if
// the visitor functions don't find a match.

using namespace Parser;

PStream::opstream&
operator<<(PStream::opstream& out, InfiniteMPO const& Op)
{
   out << Op.Operator;
   out << Op.Description;
   return out;
}

PStream::ipstream&
operator>>(PStream::ipstream& in, InfiniteMPO& Op)
{
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
      return "TriangularMPO";
   if (this->is_product())
      return "ProductMPO";
   return "unknown";
}

bool
InfiniteMPO::is_triangular() const
{
   return bool(boost::get<TriangularMPO>(&this->op()));
}

TriangularMPO const&
InfiniteMPO::as_triangular_mpo() const
{
   return boost::get<TriangularMPO>(this->op());
}
      
bool
InfiniteMPO::is_product() const
{
   return bool(boost::get<ProductMPO>(&this->op()));
}

ProductMPO const&
InfiniteMPO::as_product_mpo() const
{
   return boost::get<ProductMPO>(this->op());
}

bool
InfiniteMPO::is_complex() const
{
   return bool(boost::get<std::complex<double> >(&this->op()));
}
 
std::complex<double>
InfiniteMPO::as_complex() const
{
   return boost::get<std::complex<double> >(this->op());
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
   return boost::apply_visitor(binary_addition<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO operator-(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_subtraction<InfiniteMPOElement>(), x.op(), y.op());
}


InfiniteMPO operator-(InfiniteMPO const& x)
{
   return boost::apply_visitor(negate_element<InfiniteMPOElement>(), x.op());
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
   return boost::apply_visitor(ternary_product_q<InfiniteMPOElement>(q), x.op(), y.op());
}
   
InfiniteMPO prod(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_multiplication<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO operator*(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_multiplication<InfiniteMPOElement>(), x.op(), y.op());
}

// dot product - takes into account the multiplicity to rescale the result
InfiniteMPO dot(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_dot_product<InfiniteMPOElement>(), x.op(), y.op());
}

// inner product - equivalent to dot(adjoint(x),y)
InfiniteMPO inner(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_inner_product<InfiniteMPOElement>(), x.op(), y.op());
}

// cross product (if it exists)
InfiniteMPO cross(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_cross_product<InfiniteMPOElement>(), x.op(), y.op());
}

InfiniteMPO outer(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_outer_product<InfiniteMPOElement>(), x.op(), y.op());
}

// power of an operator.  Requires n > 1.
InfiniteMPO pow(InfiniteMPO const& x, int n)
{
   return boost::apply_visitor(unary_power<InfiniteMPOElement>(n), x.op());
}

InfiniteMPO pow(InfiniteMPO const& x, InfiniteMPO const& y)
{
   return boost::apply_visitor(binary_power<InfiniteMPOElement>(), x.op(), y.op());
}

// Exponential - only defined for complex
InfiniteMPO exp(InfiniteMPO const& x)
{
   return boost::apply_visitor(ElementExp<InfiniteMPOElement>(), x.op());
}

// Conjugate
InfiniteMPO conj(InfiniteMPO const& x)
{
   return boost::apply_visitor(ElementConj<InfiniteMPOElement>(), x.op());
}

// Adjoint
InfiniteMPO adjoint(InfiniteMPO const& x)
{
   return boost::apply_visitor(ElementAdjoint<InfiniteMPOElement>(), x.op());
}

// Inverse Adjoint
InfiniteMPO inv_adjoint(InfiniteMPO const& x)
{
   return boost::apply_visitor(ElementInvAdjoint<InfiniteMPOElement>(), x.op());
}
