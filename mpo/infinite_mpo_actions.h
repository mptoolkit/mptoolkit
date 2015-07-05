// -*- C++ -*- $Id$

// Specializations of the visitor actions for InfiniteMPO types

#if !defined(MPTOOLKIT_MPO_INFINITE_MPO_ACTIONS_H)
#define MPTOOLKIT_MPO_INFINITE_MPO_ACTIONS_H

#include "mpo/infinite_mpo.h"
#include "parser/visitor_actions.h"
#include <cmath>

namespace Parser
{

inline
std::string name_of(std::complex<double>)
{
   return "complex";
}

inline
std::string name_of(TriangularMPO const&)
{
   return "TriangularMPO";
}

inline
std::string name_of(ProductMPO const&)
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

   InfiniteMPOElement operator()(ProductMPO const&) const
   {
      PANIC("exp() is not implemented for ProductMPO's!");
      return ProductMPO();
   }

   InfiniteMPOElement operator()(TriangularMPO const&) const
   {
      PANIC("exp() is not implemented for TriangularMPO's!");
      return TriangularMPO();
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

   InfiniteMPOElement operator()(ProductMPO const&) const
   {
      PANIC("Unary negation is not defined for ProductMPO");
      return InfiniteMPOElement(0.0);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x) const
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

   InfiniteMPOElement operator()(TriangularMPO const& x) const
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

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {

      int i = int(round(y.real()));
      CHECK(std::norm(complex(double(i)) - y) < 1E-7)
	 ("cannot take fractional power of a TriangularMPO");
      CHECK(i >= 0)("cannot take negative power of a TriangularMPO");
      return i == 0 ? InfiniteMPOElement(complex(1.0,0.0)) : pow(x, i);
   }

   InfiniteMPOElement operator()(ProductMPO const& x, complex const& y) const
   {

      int i = int(round(y.real()));
      CHECK(std::norm(complex(double(i)) - y) < 1E-7)
	 ("cannot take fractional power of a ProductMPO");
      CHECK(i >= 0)("cannot take negative power of a ProductMPO");
      return i == 0 ? InfiniteMPOElement(complex(1.0,0.0)) : pow(x, i);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate exponent with an operator")(typeid(y).name());
      return InfiniteMPOElement(complex(0));
   }
};

template <>
struct binary_addition<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   { 
      return InfiniteMPOElement(x+y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(x+y);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Addition is not defined for "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
   }
};

template <>
struct binary_subtraction<InfiniteMPOElement> : boost::static_visitor<InfiniteMPOElement>
{
   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   { 
      return InfiniteMPOElement(x-y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(x-y);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Subtraction is not defined for "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   // Not defined for:
   // complex*product
   // triangular*product
   // product*complex
   // product*triangular

   InfiniteMPOElement operator()(complex const&, ProductMPO const&) const
   {
      PANIC("Cannot multiply ProductMPO by a scalar");
      return complex(0.0);
   }

   InfiniteMPOElement operator()(ProductMPO const&, complex const&) const
   {
      PANIC("Cannot multiply ProductMPO by a scalar");
      return complex(0.0);
   }

   InfiniteMPOElement operator()(TriangularMPO const&, ProductMPO const&) const
   {
      PANIC("Cannot multiply ProductMPO and a TriangularMPO");
      return complex(0.0);
   }

   InfiniteMPOElement operator()(ProductMPO const&, TriangularMPO const&) const
   {
      PANIC("Cannot multiply TriangularMPO and a ProductMPO");
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, complex const& y) const
   {
      return InfiniteMPOElement(x/y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement((1.0/y) * x);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T& x, U const& y) const
   {
      PANIC("Cannot divide "+name_of(x)+" by "+name_of(y));
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(dot(x,y));
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      return InfiniteMPOElement(dot(x,y));
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate dot product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(inner(x,y));
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      return InfiniteMPOElement(inner(x,y));
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate inner product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      return InfiniteMPOElement(x*y);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(outer(x,y));
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      return InfiniteMPOElement(outer(x,y));
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate outer product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(cross(x,y));
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate cross product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      return complex(0.0);
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return InfiniteMPOElement(x*y - y*x);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate commutator product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      QuantumNumbers::QuantumNumber q(y.GetSymmetryList(), qStr);
      CHECK_EQUAL(q, y.TransformsAs());
      return x*y;
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      QuantumNumbers::QuantumNumber q(x.GetSymmetryList(), qStr);
      CHECK_EQUAL(q, x.TransformsAs());
      return x*y;
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      QuantumNumbers::QuantumNumber q(x.GetSymmetryList(), qStr);
      return prod(x,y,q);
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      QuantumNumbers::QuantumNumber q(x.GetSymmetryList(), qStr);
      CHECK(is_scalar(q))("Product of ProductMPO's must be scalar!")(q);
      return prod(x,y);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
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

   InfiniteMPOElement operator()(complex const& x, TriangularMPO const& y) const
   {
      CHECK_EQUAL(q, y.TransformsAs());
      return x*y;
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, complex const& y) const
   {
      CHECK_EQUAL(q, x.TransformsAs());
      return x*y;
   }

   InfiniteMPOElement operator()(TriangularMPO const& x, TriangularMPO const& y) const
   {
      return prod(x,y,q);
   }

   InfiniteMPOElement operator()(ProductMPO const& x, ProductMPO const& y) const
   {
      CHECK(is_scalar(q));
      return prod(x,y);
   }

   template <typename T, typename U>
   InfiniteMPOElement operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate product between "+name_of(x)+" and "+name_of(y));
      return complex(0.0);
   }

   QuantumNumbers::QuantumNumber q;
};

} // namespace Parser

#endif
