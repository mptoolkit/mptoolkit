// -*- C++ -*- $Id$
//
// The generic part of the 'calculator' parser for reading expressions
//

#if !defined(MPTOOLKIT_PARSER_PARSER_H)
#define MPTOOLKIT_PARSER_PARSER_H

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_chset.hpp>
#include <boost/spirit/include/classic_symbols.hpp>
#include <boost/spirit/include/classic_refactoring.hpp>
#include <boost/spirit/include/classic_distinct.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <complex>
#include "common/math_const.h"
#include "quantumnumbers/quantumnumber.h"  // for the ternary_prod implementation

namespace Parser
{

using namespace boost::spirit::classic;

typedef std::complex<double> complex;

distinct_parser<> const keyword_p("0-9a-zA-Z_");

distinct_directive<> const keyword_d("0-9a-zA-Z_");

// apply a unary math function to an element
template <typename Func>
struct unary_math : boost::static_visitor<complex>
{
   unary_math() {}
   unary_math(Func f_) : f(f_) {}

   complex operator()(complex c) const
   {
      return f(c);
   }

   template <typename T>
   complex operator()(T const& t) const
   {
      PANIC("Cannot apply a mathematical function to this object")(typeid(t).name());
      return complex();
   }

   Func f;
};

template <typename element_type, typename Visitor>
struct apply_unary_math
{
   apply_unary_math(Visitor const& f_) : f(f_) {}

   typedef element_type result_type;
   typedef element_type argument_type;

   element_type operator()(element_type const& x) const
   {
      return element_type(boost::apply_visitor(f, x));
   }

   Visitor f;
};

template <typename element_type, typename F>
inline
apply_unary_math<element_type, F> make_apply_unary_math(F const& f)
{
   return apply_unary_math<element_type, F>(f);
}

template <typename element_type, typename F>
inline
apply_unary_math<element_type, unary_math<F> > make_unary_math(F const& f)
{
   return apply_unary_math<element_type, unary_math<F> >(unary_math<F>(f));
}

//
// pre-defined symbols
//

struct constants : symbols<complex>
{
   constants()
   {
      add
         ("pi", math_const::pi)
         ("e", math_const::e)
         ;
   }
};

// pre-defined unary operators

// can't find a direct way of wrapping the std:: overloaded functions :(
inline
complex csin(complex x)
{
   return sin(x);
}

inline
complex ccos(complex x)
{
   return cos(x);
}

inline
complex cexp(complex x)
{
   return exp(x);
}

inline
complex cln(complex x)
{
   return log(x);
}

inline
complex clog10(complex x)
{
   return log10(x);
}

inline
complex csqrt(complex x)
{
   return sqrt(x);
}


//typedef boost::function<element_type(element_type)> unary_func_type;

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
struct unary_funcs : symbols<boost::function<element_type(element_type)> >
{
   unary_funcs()
   {
      this->add.operator()
         ("sin", make_unary_math<element_type>(&csin))
         ("cos", make_unary_math<element_type>(&ccos))
         ("exp", make_apply_unary_math<element_type>(ElementExp<element_type>()))
         ("ln", make_unary_math<element_type>(&cln))
         ("log", make_unary_math<element_type>(&cln))
         ("log10", make_unary_math<element_type>(&clog10))
         ("sqrt", make_unary_math<element_type>(&csqrt))
         ("conj", make_apply_unary_math<element_type>(ElementConj<element_type>()))
         ("ad", make_apply_unary_math<element_type>(ElementAdjoint<element_type>()))
         ("adjoint", make_apply_unary_math<element_type>(ElementAdjoint<element_type>()))
         ("sq", make_apply_unary_math<element_type>(ElementSq<element_type>()))
         ;
   }
};

//
// binary functions
//

template <typename element_type, typename Visitor>
struct apply_binary_math
{
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;

   apply_binary_math(Visitor const& f_) : f(f_) {}

   typedef element_type result_type;
   typedef element_type first_argument_type;
   typedef element_type second_argument_type;

   element_type operator()(element_type const& x, element_type const& y) const
   {
      return element_type(boost::apply_visitor(f, x, y));
   }

   Visitor f;
};

template <typename element_type, typename F>
inline
apply_binary_math<element_type, F>
make_apply_binary_math(F const& f)
{
   return apply_binary_math<element_type, F>(f);
}

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
// cross_product: outer product of tensor operators.
// We choose among the possible transform_targets the
// quantum number with the largest degree.
template <typename element_type>
struct binary_cross_product : boost::static_visitor<element_type>
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
      return cross(x,y);
   }
};

// the standard binary functions are put in a table

template <typename element_type>
struct binary_funcs : symbols<boost::function<element_type(element_type, element_type)> >
{
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;
   binary_funcs()
   {
      this->add.operator()
         ("dot", make_apply_binary_math<element_type>(binary_dot_product<element_type>()))
         ("inner", make_apply_binary_math<element_type>(binary_inner_product<element_type>()))
         ("outer", make_apply_binary_math<element_type>(binary_outer_product<element_type>()))
	 ("cross", make_apply_binary_math<element_type>(binary_cross_product<element_type>()))
         ;
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

//
//  Semantic actions
//

struct push_identifier
{
   push_identifier(std::stack<std::string>& eval_)
      : eval(eval_) {}

   void operator()(char const* str, char const* end) const
   {
      eval.push(std::string(str, end));
   }

   std::stack<std::string>& eval;
};

template <typename element_type>
struct push_int
{
    push_int(std::stack<element_type>& eval_)
    : eval(eval_) {}

    void operator()(char const* str, char const* /*end*/) const
    {
        long n = strtol(str, 0, 10);
        eval.push(n);
        std::cout << "push\t" << long(n) << std::endl;
    }

    std::stack<element_type >& eval;
};

template <typename element_type>
struct push_real
{
    push_real(std::stack<element_type >& eval_)
    : eval(eval_) {}

   void operator()(double n) const
   {
        eval.push(n);
    }

   void operator()(element_type n) const
   {
        eval.push(n);
    }

    void operator()(char const* str, char const* /*end*/) const
    {
        double n = strtod(str, 0);
        eval.push(n);
    }

    std::stack<element_type >& eval;
};

template <typename element_type>
struct push_imag
{
    push_imag(std::stack<element_type >& eval_)
    : eval(eval_) {}

   void operator()(double n) const
   {
        eval.push(element_type(complex(0.0,n)));
    }

   void operator()(char const* beg, char const* end) const
   {
      double n = strtod(beg, 0);
      eval.push(element_type(complex(0.0,n)));
    }

    std::stack<element_type >& eval;
};

template <typename element_type>
struct push_unary
{
   typedef boost::function<element_type(element_type)> unary_func_type;
   push_unary(std::stack<unary_func_type>& s_) : s(s_) {}

   void operator()(unary_func_type f) const
   {
      s.push(f);
   }

   std::stack<unary_func_type>& s;
};

template <typename element_type>
struct eval_unary
{
   typedef boost::function<element_type(element_type)> unary_func_type;
   eval_unary(std::stack<unary_func_type>& s_, std::stack<element_type >& eval_) : s(s_), eval(eval_) {}

   void operator()(char const* beg, char const* end) const
   {
      element_type n = s.top()(eval.top()); 
      eval.pop();
      s.pop();
      eval.push(n);
   }

   std::stack<unary_func_type>& s;
   std::stack<element_type >& eval;
};

template <typename element_type>
struct push_binary
{
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;
   push_binary(std::stack<binary_func_type>& s_) : s(s_) {}

   void operator()(binary_func_type f) const
   {
      s.push(f);
   }

   std::stack<binary_func_type>& s;
};

template <typename element_type>
struct eval_binary
{
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;
   eval_binary(std::stack<binary_func_type>& s_, std::stack<element_type >& eval_) : s(s_), eval(eval_) {}

   void operator()(char const* beg, char const* end) const
   {
      element_type second = eval.top();
      eval.pop();
      element_type n = s.top()(eval.top(), second); 
      eval.pop();
      s.pop();
      eval.push(n);
   }

   std::stack<binary_func_type>& s;
   std::stack<element_type >& eval;
};

template <typename element_type, typename BinaryFunc>
struct invoke_binary
{
   invoke_binary(std::stack<element_type >& eval_) : eval(eval_) {}
   invoke_binary(std::stack<element_type >& eval_, BinaryFunc f_) : eval(eval_), f(f_) {}

   void operator()(char const*, char const*) const
   {
      element_type rhs = eval.top();
      eval.pop();
      element_type lhs = eval.top();
      eval.pop();
      eval.push(boost::apply_visitor(f, lhs, rhs));
   }

   std::stack<element_type>& eval;
   BinaryFunc f;
};

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
struct do_negate
{
    do_negate(std::stack<element_type >& eval_)
    : eval(eval_) {}

    void operator()(char const*, char const*) const
    {
        element_type lhs = eval.top();
        eval.pop();
        eval.push(boost::apply_visitor(negate_element<element_type>(), lhs));
    }

    std::stack<element_type >& eval;
};

template <typename element_type>
struct push_prod
{
   push_prod(std::stack<std::string>& identifier_stack_, std::stack<element_type >& eval_)
      : identifier_stack(identifier_stack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      element_type op2 = eval.top();
      eval.pop();
      element_type op1 = eval.top();
      eval.pop();
      std::string q = identifier_stack.top();
      identifier_stack.pop();

      eval.push(boost::apply_visitor(ternary_product<element_type>(q), op1, op2));
   }

   std::stack<std::string>& identifier_stack;
   std::stack<element_type>& eval;
};




} //namespace Parser

#endif
