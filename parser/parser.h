// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// parser/parser.h
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2024 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
//
// The generic part of the 'calculator' parser for reading expressions
//

#if !defined(MPTOOLKIT_PARSER_PARSER_H)
#define MPTOOLKIT_PARSER_PARSER_H

#include "common/math_const.h"
#include "quantumnumbers/quantumnumber.h"  // for the ternary_prod implementation
#include "linearalgebra/scalar.h"
#include "parser/visitor_actions.h"
#include "lattice/function.h"
#include "common/stringutil.h"

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_chset.hpp>
#include <boost/spirit/include/classic_symbols.hpp>
#include <boost/spirit/include/classic_refactoring.hpp>
#include <boost/spirit/include/classic_distinct.hpp>

//#include <boost/spirit/include/qi_core.hpp>

#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/get.hpp>
#include <boost/math/complex.hpp>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <complex>
#include <type_traits>
#include <utility>
#include <variant>

namespace Parser
{

// LookupFileGrid
// The filegrid(file, x, y, z) function reads values from a 3D array stored in a file.
// The LookupFileGrid function is used to lookup elements in this array.

void
LoadFileGrid(std::string Filename);

std::string
LookupFileGrid(std::string Filename, int x, int y, int z);

// returns a string of Size space characters
std::string Spaces(int Size);

// returns true if the parenthesis a b match, ie they form
// ( ) or [ ] or { }
inline
bool ParenthesesMatch(char a, char b)
{
   return (a == '(' && b == ')')
      || (a == '[' && b == ']')
      || (a == '{' && b == '}');
}

template <typename Iter>
void CheckParentheses(Iter beg, Iter end)
{
   std::stack<Iter> IterStack;
   Iter I = beg;
   while (I != end)
   {
      if (*I == '(' || *I == '[' || *I == '{')
      {
         IterStack.push(I);
      }
      else if (*I == ')' || *I == ']' || *I == '}')
      {
         if (IterStack.empty())
         {
            std::string s;
            s = s + "Unbalanced parentheses, extra '" + (*I) + "'";
            s = s + "\nWhile parsing string:\n" + std::string(beg, end);
            s = s + "\n" + Spaces(std::distance(beg, I)) + ColorPrompt("^");
            throw ParserError(s);
         }
         else
         {
            if (ParenthesesMatch(*IterStack.top(), *I))
            {
               // Found the matching parentheses
               IterStack.pop();
            }
            else
            {
               std::string s;
               s = s + "Parenthesis mismatch, opening '" + (*IterStack.top())
                  + "' closes with '" + (*I) + "'";
               s = s + "\nWhile parsing string:\n" + std::string(beg, end);
               s = s + "\n" + Spaces(std::distance(beg, IterStack.top())) + "^";
               s = s + Spaces(std::distance(IterStack.top(), I)-1) + "^";
               throw ParserError(s);
            }
         }
      }
      ++I;
   }
   if (!IterStack.empty())
   {
      std::string s;
      s = s + "Unbalanced parenthesis, '" + (*IterStack.top())
         + "' has no closing bracket";
      s = s + "\nWhile parsing string:\n" + std::string(beg, end);
      s = s + "\n" + Spaces(std::distance(beg, IterStack.top())) + "^";
      throw ParserError(s);
   }
}

inline
void CheckParentheses(char const* beg, char const* end)
{
   std::stack<char const*> IterStack;
   char const* I = beg;
   while (I != end)
   {
      if (*I == '(' || *I == '[' || *I == '{')
      {
         IterStack.push(I);
      }
      else if (*I == ')' || *I == ']' || *I == '}')
      {
         if (IterStack.empty())
         {
            throw ParserError::AtPosition(std::string("Unbalanced parenthesis, extra '")
                                          + ColorHighlight(*I) + "'",
                                          I);
         }
         else
         {
            if (ParenthesesMatch(*IterStack.top(), *I))
            {
               // Found the matching parentheses
               IterStack.pop();
            }
            else
            {
               throw ParserError::AtRange(std::string("Parenthesis mismatch, opening '")
                                          + ColorHighlight(*IterStack.top())
                                          + "' closes with '" + ColorHighlight(*I) + "'",
                                          IterStack.top(), I+1);
            }
         }
      }
      ++I;
   }
   if (!IterStack.empty())
   {
      throw ParserError::AtPosition(std::string("Unbalanced parenthesis, '")
                                    + ColorHighlight(*IterStack.top())
                                    + "' has no closing bracket",
                                    IterStack.top());
   }
}


using namespace boost::spirit::classic;

distinct_parser<> const keyword_p("0-9a-zA-Z_");

distinct_directive<> const keyword_d("0-9a-zA-Z_");

typedef std::complex<double> complex;

template <typename Visitor, typename = void>
struct has_visit_result_type : std::false_type
{};

template <typename Visitor>
struct has_visit_result_type<Visitor, std::void_t<typename std::decay_t<Visitor>::result_type>> : std::true_type
{};

template <typename Visitor, typename... Variants>
decltype(auto) apply_std_variant_visitor(Visitor&& visitor, Variants&&... variants)
{
   if constexpr (has_visit_result_type<Visitor>::value)
   {
      using result_type = typename std::decay_t<Visitor>::result_type;
      return std::visit([&visitor](auto&&... args) -> result_type
                        {
                           return std::forward<Visitor>(visitor)(std::forward<decltype(args)>(args)...);
                        },
                        std::forward<Variants>(variants)...);
   }
   else
   {
      return std::visit(std::forward<Visitor>(visitor), std::forward<Variants>(variants)...);
   }
}

template <typename Visitor, typename... T>
decltype(auto) apply_variant_visitor(Visitor&& visitor, std::variant<T...>& variant)
{
   return apply_std_variant_visitor(std::forward<Visitor>(visitor), variant);
}

template <typename Visitor, typename... T>
decltype(auto) apply_variant_visitor(Visitor&& visitor, std::variant<T...> const& variant)
{
   return apply_std_variant_visitor(std::forward<Visitor>(visitor), variant);
}

template <typename Visitor, typename... T1, typename... T2>
decltype(auto) apply_variant_visitor(Visitor&& visitor, std::variant<T1...>& x, std::variant<T2...>& y)
{
   return apply_std_variant_visitor(std::forward<Visitor>(visitor), x, y);
}

template <typename Visitor, typename... T1, typename... T2>
decltype(auto) apply_variant_visitor(Visitor&& visitor, std::variant<T1...> const& x, std::variant<T2...> const& y)
{
   return apply_std_variant_visitor(std::forward<Visitor>(visitor), x, y);
}

template <typename Visitor, typename... T1, typename... T2>
decltype(auto) apply_variant_visitor(Visitor&& visitor, std::variant<T1...>& x, std::variant<T2...> const& y)
{
   return apply_std_variant_visitor(std::forward<Visitor>(visitor), x, y);
}

template <typename Visitor, typename... T1, typename... T2>
decltype(auto) apply_variant_visitor(Visitor&& visitor, std::variant<T1...> const& x, std::variant<T2...>& y)
{
   return apply_std_variant_visitor(std::forward<Visitor>(visitor), x, y);
}

template <typename Visitor, typename Variant>
decltype(auto) apply_variant_visitor(Visitor&& visitor, Variant&& variant)
{
   return boost::apply_visitor(std::forward<Visitor>(visitor), std::forward<Variant>(variant));
}

template <typename Visitor, typename Variant1, typename Variant2>
decltype(auto) apply_variant_visitor(Visitor&& visitor, Variant1&& x, Variant2&& y)
{
   return boost::apply_visitor(std::forward<Visitor>(visitor),
                               std::forward<Variant1>(x),
                               std::forward<Variant2>(y));
}

template <typename T, typename... U>
T* variant_get(std::variant<U...>* variant)
{
   return std::get_if<T>(variant);
}

template <typename T, typename... U>
T const* variant_get(std::variant<U...> const* variant)
{
   return std::get_if<T>(variant);
}

template <typename T, typename... U>
T& variant_get(std::variant<U...>& variant)
{
   return std::get<T>(variant);
}

template <typename T, typename... U>
T const& variant_get(std::variant<U...> const& variant)
{
   return std::get<T>(variant);
}

template <typename T, typename Variant>
T* variant_get(Variant* variant)
{
   return boost::get<T>(variant);
}

template <typename T, typename Variant>
T const* variant_get(Variant const* variant)
{
   return boost::get<T>(variant);
}

template <typename T, typename Variant>
T& variant_get(Variant& variant)
{
   return boost::get<T>(variant);
}

template <typename T, typename Variant>
T const& variant_get(Variant const& variant)
{
   return boost::get<T>(variant);
}

// name_of function, for getting a friendly name of components of element_type's

template <typename T>
std::string name_of(T const&);

struct NameOf : public boost::static_visitor<std::string>
{
   template <typename T>
   std::string operator()(T const& x) const
   {
      return name_of(x);
   }
};

template <BOOST_VARIANT_ENUM_PARAMS(typename T)>
inline
std::string
name_of(boost::variant<BOOST_VARIANT_ENUM_PARAMS(T) > const& x)
{
   return apply_variant_visitor(NameOf(), x);
}

template <typename... T>
inline
std::string
name_of(std::variant<T...> const& x)
{
   return apply_variant_visitor(NameOf(), x);
}

template <>
inline
std::string name_of<complex>(std::complex<double> const&)
{
   return "complex";
}

template <typename ElementType>
int pop_int(std::stack<ElementType>& eval)
{
   if (eval.empty())
      throw ParserError("expected an integer, but the stack is empty!");
   if (!variant_get<complex>(&eval.top()))
      throw ParserError("expected an integer, got a " + name_of(eval.top()));
   complex x = variant_get<complex>(eval.top());
   eval.pop();
   return as_int(x);
}

template <typename ElementType>
double pop_real(std::stack<ElementType>& eval)
{
   if (eval.empty())
      throw ParserError("expected an integer, but the stack is empty!");
   if (!variant_get<complex>(&eval.top()))
      throw ParserError("expected an integer, got a " + name_of(eval.top()));
   complex x = variant_get<complex>(eval.top());
   eval.pop();
   return as_real(x);
}

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
      return element_type(apply_variant_visitor(f, x));
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

template <typename element_type, typename LatticeType>
struct eval_filegrid
{
   eval_filegrid(LatticeType const& Lattice_,
                 std::stack<std::string>& IdentifierStack_,
                 std::stack<int>& NumParameterStack_,
                 std::stack<element_type>& eval_,
                 Function::ArgumentList Args_ = Function::ArgumentList())
      : Lattice(Lattice_), IdentifierStack(IdentifierStack_), NumParameterStack(NumParameterStack_),
        eval(eval_),Args(Args_)  {}

   void operator()(char const* Start, char const* End) const
   {
      auto Filename = IdentifierStack.top();
      Trim(Filename);
      IdentifierStack.pop();

      int n = NumParameterStack.top();
      NumParameterStack.pop();
      if (n < 1)
         throw ParserError::AtRange(ColorHighlight("filegrid:") + " must have at least one coordinate!", Start, End);
      if (n > 3)
         throw ParserError::AtRange(ColorHighlight("filegrid:") + " cannot have more than 3 coordinates!", Start, End);

      // We represent the coordinates as the lowest numeric index, and its weight from 0..1.
      // The linear interpolation is f(x)*xw + f(x+1)*(1-xw), where x is the coordinate and xw is the weight.
      int x = 0, y = 0, z = 0;
      double xw = 1, yw = 1, zw = 1;

      double d = pop_real(eval); --n;
      x = int(d);
      xw = 1+x-d;
      if (n > 0)
      {
         y = x;
         yw = xw;
         d = pop_real(eval); --n;
         x = int(d);
         xw = 1+x-d;
         if (n > 0)
         {
            z = y;
            zw = yw;
            y = x;
            yw = xw;
            d = pop_real(eval); --n;
            x = int(d);
            xw = 1+x-d;
         }
      }
      CHECK_EQUAL(n,0);

      element_type R = ParseOperator(Lattice, LookupFileGrid(Filename, x, y, z), Args);

      // Linear interpolation in the X direction.  R = f(x)*xw + f(x+1)*(1-xw)
      if (xw < 1)
      {
         element_type X = ParseOperator(Lattice, LookupFileGrid(Filename, x+1, y, z), Args);
         element_type Xw1(1-xw);
         X = apply_variant_visitor(binary_multiplication<element_type>(), X, Xw1);
         element_type Xw(xw);
         R = apply_variant_visitor(binary_multiplication<element_type>(), R, Xw);
         R = apply_variant_visitor(binary_addition<element_type>(), R, X);
      }

      // Linear interpolation in the Y direction.  R = f(y)*yw + f(y+1)*(1-yw)
      if (yw < 1)
      {
         element_type Y = ParseOperator(Lattice, LookupFileGrid(Filename, x, y+1, z), Args);
         // Do we need a nested linear interpolation?  If so, Y = f(x,y+1)*xw + f(x+1,y+1)*(1-xw)
         if (xw < 1)
         {
            element_type XY = ParseOperator(Lattice, LookupFileGrid(Filename, x+1, y+1, z), Args);
            element_type Xw1(1-xw);
            XY = apply_variant_visitor(binary_multiplication<element_type>(), XY, Xw1);
            element_type Xw(xw);
            Y = apply_variant_visitor(binary_multiplication<element_type>(), Y, Xw);
            Y = apply_variant_visitor(binary_addition<element_type>(), Y, XY);
         }
         element_type Yw1(1-yw);
         Y = apply_variant_visitor(binary_multiplication<element_type>(), Y, Yw1);
         element_type Yw(yw);
         R = apply_variant_visitor(binary_multiplication<element_type>(), R, Yw);
         R = apply_variant_visitor(binary_addition<element_type>(), R, Y);
      }

      // Linear interpolation in the Z direction.  R = f(z)*zw + f(z+1)*(1-zw)
      if (zw < 1)
      {
         element_type Z = ParseOperator(Lattice, LookupFileGrid(Filename, x, y, z+1), Args);
         if (yw < 1)
         {
            // nested interpolation in Y,Z directions
            element_type YZ = ParseOperator(Lattice, LookupFileGrid(Filename, x, y+1, z+1), Args);
            if (xw < 1)
            {
               // nested interpolation in X,Y,Z directions
               element_type XZ = ParseOperator(Lattice, LookupFileGrid(Filename, x+1, y, z+1), Args);
               element_type Xw1(1-xw);
               XZ = apply_variant_visitor(binary_multiplication<element_type>(), XZ, Xw1);
               element_type Xw(xw);
               Z = apply_variant_visitor(binary_multiplication<element_type>(), Z, Xw);
               Z = apply_variant_visitor(binary_addition<element_type>(), Z, XZ);

               element_type XYZ = ParseOperator(Lattice, LookupFileGrid(Filename, x+1, y+1, z+1), Args);
               XYZ = apply_variant_visitor(binary_multiplication<element_type>(), XYZ, Xw1);
               YZ = apply_variant_visitor(binary_multiplication<element_type>(), YZ, Xw);
               YZ = apply_variant_visitor(binary_addition<element_type>(), YZ, XYZ);
            }
            element_type Yw1(1-yw);
            YZ = apply_variant_visitor(binary_multiplication<element_type>(), YZ, Yw1);
            element_type Yw(yw);
            Z = apply_variant_visitor(binary_multiplication<element_type>(), Z, Yw);
            Z = apply_variant_visitor(binary_addition<element_type>(), Z, YZ);
         }
         else if (xw < 1)
         {
            // nested interpolation in X,Z directions
            element_type XZ = ParseOperator(Lattice, LookupFileGrid(Filename, x+1, y, z+1), Args);
            element_type Xw1(1-xw);
            XZ = apply_variant_visitor(binary_multiplication<element_type>(), XZ, Xw1);
            element_type Xw(xw);
            Z = apply_variant_visitor(binary_multiplication<element_type>(), Z, Xw);
            Z = apply_variant_visitor(binary_addition<element_type>(), Z, XZ);
         }
         element_type Zw1(1-zw);
         Z = apply_variant_visitor(binary_multiplication<element_type>(), Z, Zw1);
         element_type Zw(zw);
         R = apply_variant_visitor(binary_multiplication<element_type>(), R, Zw);
         R = apply_variant_visitor(binary_addition<element_type>(), R, Z);
      }

      eval.push(R);
   }

   LatticeType const& Lattice;
   std::stack<std::string>& IdentifierStack;
   std::stack<int>& NumParameterStack;
   std::stack<element_type>& eval;
   Function::ArgumentList Args;
};

//
// pre-defined symbols
//

struct constants : symbols<complex>
{
   constants()
   {
      add
         ("inf", std::numeric_limits<double>::infinity())
         ("nan", std::numeric_limits<double>::quiet_NaN())
         ("pi", math_const::pi)
         ("e", math_const::e)
         ("i", std::complex<double>(0.0, 1.0))
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
complex ctan(complex x)
{
   return tan(x);
}

inline
complex casin(complex x)
{
   return std::asin(x);
}

inline
complex cacos(complex x)
{
   return std::acos(x);
}

inline
complex catan(complex x)
{
   return std::atan(x);
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


//typedef std::function<element_type(element_type)> unary_func_type;

template <typename element_type>
struct unary_funcs : symbols<std::function<element_type(element_type)> >
{
   unary_funcs()
   {
      this->add.operator()
         ("sin", make_unary_math<element_type>(&csin))
         ("cos", make_unary_math<element_type>(&ccos))
         ("tan", make_unary_math<element_type>(&ctan))
         ("asin", make_unary_math<element_type>(&casin))
         ("acos", make_unary_math<element_type>(&cacos))
         ("atan", make_unary_math<element_type>(&catan))
         ("exp", make_apply_unary_math<element_type>(ElementExp<element_type>()))
         ("abs", make_apply_unary_math<element_type>(ElementAbs<element_type>()))
         ("ln", make_unary_math<element_type>(&cln))
         ("log", make_unary_math<element_type>(&cln))
         ("log10", make_unary_math<element_type>(&clog10))
         ("sqrt", make_unary_math<element_type>(&csqrt))
         ("conj", make_apply_unary_math<element_type>(ElementConj<element_type>()))
         ("ad", make_apply_unary_math<element_type>(ElementAdjoint<element_type>()))
         ("adjoint", make_apply_unary_math<element_type>(ElementAdjoint<element_type>()))
         ("inv_adjoint", make_apply_unary_math<element_type>(ElementInvAdjoint<element_type>()))
         ("sq", make_apply_unary_math<element_type>(ElementSq<element_type>()))
         ;
   }
};

template <typename element_type>
struct unary_funcs_mpo : unary_funcs<element_type>
{
   unary_funcs_mpo()
   {
      this->add.operator()
      ("gauge_flip", make_apply_unary_math<element_type>(element_gauge_flip<element_type>()))
      ;
   }
};
//
// binary functions
//

template <typename element_type, typename Visitor>
struct apply_binary_math
{
   typedef std::function<element_type(element_type, element_type)> binary_func_type;

   apply_binary_math(Visitor const& f_) : f(f_) {}

   typedef element_type result_type;
   typedef element_type first_argument_type;
   typedef element_type second_argument_type;

   element_type operator()(element_type const& x, element_type const& y) const
   {
      return element_type(apply_variant_visitor(f, x, y));
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

// the standard binary functions are put in a table

template <typename element_type>
struct binary_funcs : symbols<std::function<element_type(element_type, element_type)> >
{
   typedef std::function<element_type(element_type, element_type)> binary_func_type;
   binary_funcs()
   {
      this->add.operator()
         ("pow", make_apply_binary_math<element_type>(binary_power<element_type>()))
         ("fmod", make_apply_binary_math<element_type>(binary_fmod<element_type>()))
         ("dot", make_apply_binary_math<element_type>(binary_dot_product<element_type>()))
         ("inner", make_apply_binary_math<element_type>(binary_inner_product<element_type>()))
         ("outer", make_apply_binary_math<element_type>(binary_outer_product<element_type>()))
         ("cross", make_apply_binary_math<element_type>(binary_cross_product<element_type>()))
         ;
   }
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
        //std::cout << "push\t" << long(n) << std::endl;
    }

    std::stack<element_type >& eval;
};

template <typename element_type>
struct push_value
{
   push_value(std::stack<element_type >& eval_)
    : eval(eval_) {}

   void operator()(double n) const
   {
      element_type c(n);
      eval.push(element_type(n));
   }

   void operator()(std::complex<double> n) const
   {
      eval.push(n);
   }

   template <typename T>
   void operator()(T const& x) const
   {
      eval.push(x);
   }

   void operator()(char const* str, char const* /*end*/) const
   {
      double n = strtod(str, 0);
      eval.push(n);
   }

   std::stack<element_type>& eval;
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
   typedef std::function<element_type(element_type)> unary_func_type;
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
   typedef std::function<element_type(element_type)> unary_func_type;
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
   typedef std::function<element_type(element_type, element_type)> binary_func_type;
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
   typedef std::function<element_type(element_type, element_type)> binary_func_type;
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

   void operator()(char const* a, char const* b) const
   {
      try
      {
         element_type rhs = eval.top();
         eval.pop();
         element_type lhs = eval.top();
         eval.pop();
         eval.push(apply_variant_visitor(f, lhs, rhs));
      }
      catch (ParserError const& p)
      {
         throw ParserError::AtPosition(p, a);
      }
      catch (std::exception const& p)
      {
         throw ParserError::AtPosition(p, a);
      }
      catch (...)
      {
         throw;
      }
   }

   std::stack<element_type>& eval;
   BinaryFunc f;
};

template <typename element_type>
struct do_negate
{
   do_negate(std::stack<element_type >& eval_)
       : eval(eval_) {}

   void operator()(char const* a, char const*) const
   {
      try
      {
         element_type lhs = eval.top();
         eval.pop();
         eval.push(apply_variant_visitor(negate_element<element_type>(), lhs));
      }
      catch (ParserError const& p)
      {
         throw ParserError::AtPosition(p, a);
      }
      catch (std::exception const& p)
      {
         throw ParserError::AtPosition(p, a);
      }
      catch (...)
      {
         throw;
      }
   }

   std::stack<element_type >& eval;
};

template <typename element_type>
struct push_prod
{
   push_prod(std::stack<std::string>& identifier_stack_, std::stack<element_type >& eval_)
      : identifier_stack(identifier_stack_), eval(eval_) {}

   void operator()(char const* a, char const*) const
   {
      try
      {
         element_type op2 = eval.top();
         eval.pop();
         element_type op1 = eval.top();
         eval.pop();
         std::string q = identifier_stack.top();
         identifier_stack.pop();

         eval.push(apply_variant_visitor(ternary_product<element_type>(q), op1, op2));
      }
      catch (ParserError const& p)
      {
         throw ParserError::AtPosition(p, a);
      }
      catch (std::exception const& p)
      {
         throw ParserError::AtPosition(p, a);
      }
      catch (...)
      {
         throw;
      }
   }

   std::stack<std::string>& identifier_stack;
   std::stack<element_type>& eval;
};

// functions

// Adds a function name to the function stack, and adds a new (initially empty) ParameterList
// to the parameter stack
struct push_function
{   push_function(std::stack<std::string>& FuncStack_, std::stack<Function::ParameterList>& ParamStack_)
      : FuncStack(FuncStack_), ParamStack(ParamStack_) {}

   void operator()(char const* str, char const* end) const
   {
      FuncStack.push(std::string(str, end));
      ParamStack.push(Function::ParameterList());
   }

   std::stack<std::string>& FuncStack;
   std::stack<Function::ParameterList>& ParamStack;
};

template <typename element_type>
struct push_parameter
{
   push_parameter(std::stack<element_type>& eval_,
                  std::stack<Function::ParameterList>& ParamStack_)
      : eval(eval_), ParamStack(ParamStack_) {}

   void operator()(char const*, char const*) const
   {
      ParamStack.top().push_back(Function::Parameter(variant_get<complex>(eval.top())));
      eval.pop();
   }

   std::stack<element_type>& eval;
   std::stack<Function::ParameterList>& ParamStack;
};

template <typename element_type>
struct push_named_parameter
{
   push_named_parameter(std::stack<element_type>& eval_,
                        std::stack<std::string>& IdentifierStack_,
                        std::stack<Function::ParameterList>& ParamStack_)
      : eval(eval_), IdentifierStack(IdentifierStack_), ParamStack(ParamStack_) {}

   void operator()(char const*, char const*) const
   {
      ParamStack.top().push_back(Function::Parameter(IdentifierStack.top(),
                                                     variant_get<complex>(eval.top())));
      IdentifierStack.pop();
      eval.pop();
   }

   std::stack<element_type>& eval;
   std::stack<std::string>& IdentifierStack;
   std::stack<Function::ParameterList>& ParamStack;
};

} //namespace Parser

#endif
