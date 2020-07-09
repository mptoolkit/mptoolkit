// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mpo/infinite_mpo_actions.h
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

// Specializations of the visitor actions for InfiniteMPO types
// TODO: The exception here shouldn't really be a ParserError, it is more like
// a type mismatch error.

#if !defined(MPTOOLKIT_MPO_INFINITE_MPO_ACTIONS_H)
#define MPTOOLKIT_MPO_INFINITE_MPO_ACTIONS_H

#include "parser/parser.h"
#include "mpo/infinite_mpo.h"
#include "parser/visitor_actions.h"
#include <cmath>

namespace Parser
{

template <>
inline
std::string name_of<ZeroMPO>(ZeroMPO const&)
{
   return "ZeroMPO";
}

template <>
inline
std::string name_of<BasicTriangularMPO>(BasicTriangularMPO const&)
{
   return "BasicTriangularMPO";
}

template <>
inline
std::string name_of<ProductMPO>(ProductMPO const&)
{
   return "ProductMPO";
}

template <>
struct ElementExp<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex c) const
   {
      return std::exp(c);
   }

   InfiniteMPOElement operator()(ZeroMPO const&) const
   {
      return std::complex<double>(1.0);
   }

   InfiniteMPOElement operator()(ProductMPO const&) const
   {
      throw ParserError("exp() is not implemented for ProductMPO's!");
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const&) const
   {
      throw ParserError("exp() is not implemented for BasicTriangularMPO's!");
   }
};

template <>
struct ElementAbs<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex c) const
   {
      return std::norm(c);
   }

   InfiniteMPOElement operator()(ZeroMPO const& c) const
   {
      return c;
   }

   InfiniteMPOElement operator()(ProductMPO const&) const
   {
      throw ParserError("abs() is not implemented for ProductMPO's!");
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const&) const
   {
      throw ParserError("abs() is not implemented for BasicTriangularMPO's!");
   }
};

// Default for ElementSq is OK
// Default for ElementAdjoint is OK
// Default for ElementInvAdjoint is OK
// Default for ElementConj is OK

template <>
struct negate_element<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex c) const
   {
      return -c;
   }

   InfiniteMPOElement operator()(ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   InfiniteMPOElement operator()(ProductMPO const&) const
   {
      throw ParserError("Unary negation is not defined for ProductMPO");
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x) const
   {
      return -x;
   }
};

#if 0
// not needed - the default version is OK
template <>
struct unary_power<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   unary_power(int n_) : n(n_) {}

   InfiniteMPOElement operator()(complex c) const
   {
      return std::pow(c, n);
   }

   InfiniteMPOElement operator()(ProductMPO const& x) const
   {
      return pow(x, n);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x) const
   {
      return pow(x, n);
   }

   int n;
};
#endif

template <>
struct binary_power<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return std::pow(x,y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {

      int i = int(round(y.real()));
      if (std::norm(complex(double(i)) - y) > 1E-7)
         throw ParserError("cannot take a fractional or complex power " + formatting::format_complex(y)
                           + " of a BasicTriangularMPO");
      if (i < 0)
         throw ParserError("cannot take negative power " + boost::lexical_cast<std::string>(i)
                           + " of a BasicTriangularMPO");
      return i == 0 ? InfiniteMPOElement(complex(1.0,0.0)) : pow(x, i);
   }

   InfiniteMPOElement operator()(ProductMPO const& x, complex const& y) const
   {

      int i = int(round(y.real()));
      if (std::norm(complex(double(i)) - y) > 1E-7)
         throw ParserError("cannot take a fractional or complex power " + formatting::format_complex(y)
                           + " of a ProductMPO");
      if (i < 0)
         throw ParserError("cannot take negative power " + boost::lexical_cast<std::string>(i)
                           + " of a ProductMPO");
      return i == 0 ? InfiniteMPOElement(complex(1.0,0.0)) : pow(x, i);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("pow(x,y) is not defined for operator y");
   }
};

template <>
struct binary_addition<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x+y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(x+y);
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const& y) const
   {
      return y;
   }

   template <typename T>
   InfiniteMPOElement operator()(T const& x, ZeroMPO const&) const
   {
      return x;
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Addition is not defined for "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct binary_subtraction<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x-y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(x-y);
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   InfiniteMPOElement operator()(ZeroMPO const&, BasicTriangularMPO const& y) const
   {
      return -y;
   }

   InfiniteMPOElement operator()(ZeroMPO const&, complex const& y) const
   {
      return -y;
   }

   template <typename T>
   InfiniteMPOElement operator()(T const& x, ZeroMPO const&) const
   {
      return x;
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Subtraction is not defined for "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct binary_multiplication<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // Multiplication is defined for:
   // complex*complex
   // complex*triangular
   // triangular*complex
   // triangular*triangular
   // product*product
   // anything*zero
   // zero*anything

   // Not defined for:
   // complex*product
   // triangular*product
   // product*complex
   // product*triangular

   InfiniteMPOElement operator()(complex const&, ProductMPO const&) const
   {
      throw ParserError("Cannot multiply ProductMPO by a scalar");
   }

   InfiniteMPOElement operator()(ProductMPO const&, complex const&) const
   {
      throw ParserError("Cannot multiply ProductMPO by a scalar");
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const&, ProductMPO const&) const
   {
      throw ParserError("Cannot multiply ProductMPO and BasicTriangularMPO");
   }

   InfiniteMPOElement operator()(ProductMPO const&, BasicTriangularMPO const&) const
   {
      throw ParserError("Cannot multiply ProductMPO and BasicTriangularMPO");
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      return InfiniteMPOElement(x*y);
   }
};

template <>
struct binary_division<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // division is defined for:
   // complex/complex
   // triangular/complex
   // zero/anything

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x/y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement((1.0/y) * x);
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot divide "+name_of(x)+" by "+name_of(y));
      return complex(0.0);
   }
};

template <>
struct binary_dot_product<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // dot() is defined for:
   // complex,complex
   // triangular,complex
   // complex,triangular
   // triangular,triangular
   // product,product

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(dot(x,y));
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      return InfiniteMPOElement(dot(x,y));
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot evaluate dot product between "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct binary_inner_product<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // inner() is defined for:
   // complex,complex
   // triangular,complex
   // complex,triangular
   // triangular,triangular
   // product,product

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(inner(x,y));
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      return InfiniteMPOElement(inner(x,y));
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot evaluate inner product between "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct binary_outer_product<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // outer() is defined for:
   // complex,complex
   // triangular,complex
   // complex,triangular
   // triangular,triangular
   // product,product

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(outer(x,y));
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      return InfiniteMPOElement(outer(x,y));
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot evaluate outer product between "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct binary_cross_product<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // cross() is defined for:
   // complex,complex        (result is zero)
   // triangular,complex     (result is zero)
   // complex,triangular     (result is zero)
   // triangular,triangular

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(cross(x,y));
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot evaluate cross product between "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct binary_commutator<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // commutator is defined for:
   // complex,complex        (result is zero)
   // triangular,complex     (result is zero)
   // complex,triangular     (result is zero)
   // triangular,triangular

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y - y*x);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot evaluate commutator product between "+name_of(x)+" and "+name_of(y));
   }
};

template <>
struct ternary_product<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // prod(x,y,q) is defined for:
   // complex*complex        (we ignore the quantum number, since we can't parse it)
   // complex*triangular
   // triangular*complex
   // triangular*triangular
   // product*product        (only if q is scalar)

   // Not defined for:
   // triangular*product
   // product*triangular
   // complex*product
   // product*complex

   ternary_product(std::string const& qStr_) : qStr(qStr_) {}

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      QuantumNumbers::QuantumNumber q(y.GetSymmetryList(), qStr);
      if (q != y.TransformsAs())
      {
         throw ParserError("In prod(x,y,q): quantum number is not possible."
                           "\nq = " + q.ToString()
                           + "\nthe only possible quantum number in this context is "
                           + y.TransformsAs().ToString());
      }
      return x*y;
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      QuantumNumbers::QuantumNumber q(x.GetSymmetryList(), qStr);
      if (q != x.TransformsAs())
      {
         throw ParserError("In prod(x,y,q): quantum number is not possible."
                           "\nq = " + q.ToString()
                           + "\nthe only possible quantum number in this context is "
                           + x.TransformsAs().ToString());
      }
      return x*y;
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      QuantumNumbers::QuantumNumber q(x.GetSymmetryList(), qStr);
      return prod(x,y,q);
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      QuantumNumbers::QuantumNumber q(x.GetSymmetryList(), qStr);
      if (!is_scalar(q))
      {
         throw ParserError("In operator prod(x,y,q): the quantum number must be scalar."
                           "\nq = " + q.ToString());
      }
      return x*y;
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("In operator prod(x,y,q): "
                        "cannot evaluate product between "+name_of(x)+" and "+name_of(y));
   }

   std::string qStr;
};

template <>
struct ternary_product_q<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   // prod(x,y,q) is defined for:
   // complex*complex        (we ignore the quantum number, since we can't parse it)
   // complex*triangular
   // triangular*complex
   // triangular*triangular
   // product*product        (only if q is scalar)

   // Not defined for:
   // triangular*product
   // product*triangular
   // complex*product
   // product*complex

   ternary_product_q(QuantumNumbers::QuantumNumber const& q_) : q(q_) {}

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      CHECK(is_scalar(q));
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(complex const& x, BasicTriangularMPO const& y) const
   {
      CHECK_EQUAL(q, y.TransformsAs());
      return x*y;
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, complex const& y) const
   {
      CHECK_EQUAL(q, x.TransformsAs());
      return x*y;
   }

   InfiniteMPOElement operator()(BasicTriangularMPO const& x, BasicTriangularMPO const& y) const
   {
      return prod(x,y,q);
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      CHECK(is_scalar(q));
      return prod(x,y);
   }

   InfiniteMPOElement operator()(ZeroMPO const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(ZeroMPO const&, T const&) const
   {
      return ZeroMPO();
   }

   template <typename T>
   InfiniteMPOElement operator()(T const&, ZeroMPO const&) const
   {
      return ZeroMPO();
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      throw ParserError("Cannot evaluate product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
   }

   QuantumNumbers::QuantumNumber q;
};

} // namespace Parser

#endif
