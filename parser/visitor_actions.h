// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// parser/visitor_actions.h
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

// templated visitor actions for various unary and binary operators

#if !defined(MPTOOLKIT_PARSER_VISITOR_ACTIONS_H)
#define MPTOOLKIT_PARSER_VISITOR_ACTIONS_H

#include "quantumnumbers/quantumnumber.h"
#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <complex>
#include <boost/math/special_functions/round.hpp>
#include "common/numerics.h"
#include "common/formatting.h"

namespace Parser
{

using formatting::format_complex;

enum ShowColors { ColorNever, ColorAuto, ColorAlways };

// set the color for output, never, auto, always
// auto: show color iff std::cout is a terminal
void SetShowColors(ShowColors c);

std::string ColorHighlight(std::string s);

std::string ColorHighlight(char c);

std::string ColorError(std::string s);

std::string ColorWarning(std::string s);

std::string ColorNote(std::string s);

std::string ColorHint(std::string s);

std::string ColorPrompt(std::string s);

std::string ColorQuote(std::string s);

typedef std::complex<double> complex;

class ParserError : public std::exception
{
   public:
      explicit ParserError(std::string const& Why);

      ParserError(ParserError const& Prev, std::string const& Why);

      ParserError(std::exception const& Prev, std::string const& Why);

      // Adds a 'hint' message to the error
      void AddHint(std::string const& s) { Hint = s; }

      std::string const& hint() const { return Hint; }

      ~ParserError() noexcept { }

      virtual const char* what() const noexcept { return Msg.c_str(); }

      // named constructors

      // Add an iterator position at the point where the error occurs
      static ParserError AtPosition(std::string const& Why, char const* Position);
      static ParserError AtPosition(ParserError const& Prev, char const* Position);
      static ParserError AtPosition(std::exception const& Prev, char const* Position);

      static ParserError AtRange(std::string const& Why, char const* Start, char const* End);
      static ParserError AtRange(ParserError const& Prev, char const* Start, char const* End);
      static ParserError AtRange(std::exception const& Prev, char const* Start, char const* End);

      // finalize once we have the complete string
      static ParserError Finalize(ParserError const& Prev, std::string const& Why,
                                  char const* beg, char const* end);

      static ParserError Finalize(std::exception const& Prev, std::string const& Why,
                                  char const* beg, char const* end);

   private:
      ParserError(std::list<std::string> const& CallStack_, char const* Position,
                  std::string const& Hint_ = "");

      ParserError(std::list<std::string> const& CallStack_, char const* Position, char const* End_,
                  std::string const& Hint_ = "");
      ParserError(std::list<std::string> const& CallStack_,
                  std::string const& Why, char const* Position, char const* End_,
                  char const* beg, char const* end,
                  std::string const& Hint_ = "");

      void AssembleMessage();

      std::list<std::string> CallStack;
      std::string Hint;
      std::string Msg;

      char const* Pos;
      char const* End;
};

inline
int as_int(complex x)
{
   int j = boost::math::iround(x.real());
   if (LinearAlgebra::norm_frob(x - double(j)) > 1E-7)
       throw ParserError("expected an integer, got a real/complex number: " + format_complex(x));
   return j;
}

inline
double as_real(complex x)
{
   if (LinearAlgebra::norm_frob(x.imag()) > 1E-7)
       throw ParserError("expected a real number, got a complex number: " + format_complex(x));
   return x.real();
}

//
// unary functions
//

template <typename element_type>
struct negate_element : boost::static_visitor<element_type>
{
   template <typename T>
   element_type operator()(T const& x) const
   {
      return element_type(-x);
   }
};

template <typename element_type>
struct ElementConj : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return LinearAlgebra::conj(c);
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return conj(x);
   }
};

template <typename element_type>
struct ElementAdjoint : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return LinearAlgebra::conj(c);
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return adjoint(x);
   }
};

template <typename element_type>
struct ElementInvAdjoint : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return LinearAlgebra::conj(c);
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return inv_adjoint(x);
   }
};

template <typename element_type>
struct ElementSq : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return complex(c.real()*c.real() + c.imag()*c.imag());
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return dot(adjoint(x),x);
   }
};

template <typename element_type>
struct ElementExp : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return std::exp(c);
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return exp(x);
   }
};

template <typename element_type>
struct ElementAbs : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return std::abs(c);
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      using LinearAlgebra::abs;
      return abs(x);
   }
};

template <typename element_type>
struct element_flip_conj : boost::static_visitor<element_type>
{
   element_flip_conj() {}

   std::complex<double> operator()(std::complex<double> const& x) const
   {
      return x;
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return flip_conj(x);
   }
};

template <typename element_type, typename BasisType>
struct element_flip_conj_basis : boost::static_visitor<element_type>
{
   element_flip_conj_basis(BasisType const& r_) : r(r_) {}

   std::complex<double> operator()(std::complex<double> const& x) const
   {
      return x;
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return flip_conj(x, r);
   }

   BasisType const& r;
};

template <typename element_type>
struct unary_power : boost::static_visitor<element_type>
{
   unary_power(int n_) : n(n_) {}

   element_type operator()(complex const& x) const
   {
      return std::pow(x,n);
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      return pow(x, n);
   }

   int n;
};

//
// binary functions
//

template <typename element_type>
struct binary_addition : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x+y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(x*MakeIdentityFrom(y) + y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(x + y*MakeIdentityFrom(x));
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      return element_type(x+y);
   }
};

template <typename element_type>
struct binary_subtraction : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x-y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(x*MakeIdentityFrom(y) - y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(x - y*MakeIdentityFrom(x));
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      return element_type(x-y);
   }
};

template <typename element_type>
struct binary_multiplication : boost::static_visitor<element_type>
{
   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return element_type(x*y);
   }
};

template <typename element_type>
struct binary_division : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x/y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      PANIC("Cannot use operator as divisor")(typeid(y).name());
      return element_type(complex());
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type((1.0 / y) * x);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Cannot divide these types")(typeid(x).name())(typeid(y).name());
      return element_type(complex());
   }
};

template <typename element_type>
struct binary_modulus : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      if (y.real() == 0)
         throw ParserError("Division by zero in modulus operator");

      // make sure that the divisor is real
      if (std::abs(y.imag() * y.real()) > 1000*std::numeric_limits<double>::epsilon())
         throw ParserError("Divisor must be real for modulus operator");

      int xReal = int(std::rint(x.real()));
      int xImag = int(std::rint(x.imag()));
      int yReal = int(std::rint(y.real()));

      return std::complex<double>(numerics::divd(xReal, yReal).rem,
                                  numerics::divd(xImag, yReal).rem);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      throw ParserError("Operands of % must be numeric types");
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      throw ParserError("Operands of % must be numeric types");
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      throw ParserError("Operands of % must be numeric types");
   }
};

template <typename element_type>
struct binary_fmod : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      if (y.real() == 0)
         throw ParserError("Division by zero in modulus operator");

      // make sure that the divisor is real
      if (std::abs(y.imag() * y.real()) > 1000*std::numeric_limits<double>::epsilon())
         throw ParserError("Divisor must be real for modulus operator");

      return std::complex<double>(std::fmod(x.real(), y.real()), std::fmod(x.imag(), y.real()));
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      throw ParserError("Operands of fmod() must be numeric types");
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      throw ParserError("Operands of fmod() must be numeric types");
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      throw ParserError("Operands of fmod() must be numeric types");
   }
};

template <typename element_type>
struct binary_power : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return std::pow(x,y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      int i = as_int(y);
      if (i < 0)
         throw ParserError("cannot take negative power of an operator");
      return i == 0 ? element_type(complex(1.0,0.0)) : pow(x, i);
   }

   // a^X == exp(ln(a)*X), so we can handle the case where a is a complex number
   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return exp(std::log(x)*y);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate exponent with an operator")(typeid(y).name());
      return element_type(complex(0));
   }
};

template <typename element_type>
struct binary_commutator : boost::static_visitor<element_type>
{
   // commutator is zero whenever we have something that commutes
   element_type operator()(complex const&, complex const&) const
   {
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(T const&, complex const&) const
   {
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(complex const&, T const&) const
   {
      return complex(0.0);
   }

   // in the generic case, calculate the actual commutator
   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return commutator(x,y);
   }
};

template <typename element_type>
struct binary_dot_product : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   // assume here we have operators
   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return dot(x,y);
   }
};

template <typename element_type>
struct binary_inner_product : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(LinearAlgebra::conj(x)*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(LinearAlgebra::conj(x)*y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(adjoint(x)*y);
   }

   // assume here we have operators
   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return dot(adjoint(x),y);
   }
};

//
// outer_product: outer product of tensor operators.
// We choose among the possible transform_targets the
// quantum number with the largest degree.
template <typename element_type>
struct binary_outer_product : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return outer(x,y);
   }
};

//
// cross_product
template <typename element_type>
struct binary_cross_product : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(0.0);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(0.0);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(0.0);
   }

   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return cross(x,y);
   }
};

//
// Ternary functions
//

template <typename element_type>
struct ternary_product : boost::static_visitor<element_type>
{
   ternary_product(std::string const& q_)
      : q(q_) {}

   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      // legitimate use of TransformsAs() on MPOs
#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
      QuantumNumbers::QuantumNumber qn(x.GetSymmetryList(), q);
      return element_type(prod(x,y,qn));
#else
      return x;
#endif
   }

   std::string q;
};

template <typename element_type>
struct ternary_product_q : boost::static_visitor<element_type>
{
   ternary_product_q(QuantumNumbers::QuantumNumber const& q_)
      : q(q_) {}

   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      // legitimate use of TransformsAs() on MPOs
#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
      return element_type(prod(x,y,q));
#else
      return x;
#endif
   }

   QuantumNumbers::QuantumNumber q;
};

} // namespace Parser

#endif
