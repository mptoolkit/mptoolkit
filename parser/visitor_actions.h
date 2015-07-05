// -*- C++ -*- $Id$

// templated visitor actions for various unary and binary operators

#if !defined(MPTOOLKIT_PARSER_VISITOR_ACTIONS_H)
#define MPTOOLKIT_PARSER_VISITOR_ACTIONS_H

#include "quantumnumbers/quantumnumber.h"
#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <complex>

namespace Parser
{

typedef std::complex<double> complex;

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
struct binary_power : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return std::pow(x,y);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      int i = int(y.real()+0.5);
      CHECK(std::norm(complex(double(i)) - y) < 1E-7)("cannot take fractional power of an operator");
      CHECK(i >= 0)("cannot take negative power of an operator");
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
      return element_type(x*y - y*x);
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

