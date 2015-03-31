// -*- C++ -*- $Id$

// *** OBSOLETE ***

#include "operator-parser.h"
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
#include "matrixproduct/mpoperator.h"

////////////////////////////////////////////////////////////////////////////

using namespace boost::spirit;

typedef std::complex<double> complex;

distinct_parser<> const keyword_p("0-9a-zA-Z_");

distinct_directive<> const keyword_d("0-9a-zA-Z_");

// our element type for the calculator stack
typedef boost::variant<complex, MPOperator> element_type;

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

   complex operator()(MPOperator const& x) const
   {
      PANIC("Cannot apply a mathematical function to an MPOperator");
      return complex();
   }

   Func f;
};

template <typename Visitor>
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

template <typename F>
inline
apply_unary_math<F> make_apply_unary_math(F const& f)
{
   return apply_unary_math<F>(f);
}

template <typename F>
inline
apply_unary_math<unary_math<F> > make_unary_math(F const& f)
{
   return apply_unary_math<unary_math<F> >(unary_math<F>(f));
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
complex csin(complex x)
{
   return sin(x);
}

complex ccos(complex x)
{
   return cos(x);
}

complex cexp(complex x)
{
   return exp(x);
}

complex cln(complex x)
{
   return log(x);
}

complex clog10(complex x)
{
   return log10(x);
}

complex csqrt(complex x)
{
   return sqrt(x);
}


typedef boost::function<element_type(element_type)> unary_func_type;

struct ElementConj : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return LinearAlgebra::conj(c);
   }

   element_type operator()(MPOperator const& x) const
   {
      return conj(x);
   }
};

struct ElementAdjoint : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return LinearAlgebra::conj(c);
   }

   element_type operator()(MPOperator const& x) const
   {
      return adjoint(x);
   }
};

struct unary_funcs : symbols<unary_func_type>
{
   unary_funcs()
   {
      add
         ("sin", make_unary_math(&csin))
         ("cos", make_unary_math(&ccos))
         ("exp", make_unary_math(&cexp))
         ("ln", make_unary_math(&cln))
         ("log", make_unary_math(&cln))
         ("log10", make_unary_math(&clog10))
         ("sqrt", make_unary_math(&csqrt))
         ("conj", make_apply_unary_math(ElementConj()))
         ("ad", make_apply_unary_math(ElementAdjoint()))
         ("adjoint", make_apply_unary_math(ElementAdjoint()))
         ;
   }
};

//
// binary functions
//

typedef boost::function<element_type(element_type, element_type)> binary_func_type;

template <typename Visitor>
struct apply_binary_math
{
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

template <typename F>
inline
apply_binary_math<F> make_apply_binary_math(F const& f)
{
   return apply_binary_math<F>(f);
}

struct binary_addition : boost::static_visitor<element_type>
{
   binary_addition(OperatorList const& OpList_) : OpList(OpList_) {}

   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x+y);
   }

   element_type operator()(MPOperator const& x, MPOperator const& y) const
   {
      return element_type(x+y);
   }

   element_type operator()(complex const& x, MPOperator const& y) const
   {
      return element_type(x*OpList["I"] + y);
   }

   element_type operator()(MPOperator const& x, complex const& y) const
   {
      return element_type(x + y*OpList["I"]);
   }

   OperatorList const& OpList;
};

struct binary_subtraction : boost::static_visitor<element_type>
{
   binary_subtraction(OperatorList const& OpList_) : OpList(OpList_) {}

   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x-y);
   }

   element_type operator()(MPOperator const& x, MPOperator const& y) const
   {
      return element_type(x-y);
   }

   element_type operator()(complex const& x, MPOperator const& y) const
   {
      return element_type(x*OpList["I"] - y);
   }

   element_type operator()(MPOperator const& x, complex const& y) const
   {
      return element_type(x - y*OpList["I"]);
   }

   OperatorList const& OpList;
};

struct binary_multiplication : boost::static_visitor<element_type>
{
   template <typename T1, typename T2>
   element_type operator()(T1& x, T2& y) const
   {
      return element_type(x*y);
   }
};

struct binary_division : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x/y);
   }

   element_type operator()(MPOperator const& x, MPOperator const& y) const
   {
      return element_type(x-y);
   }

   element_type operator()(complex const& x, MPOperator const& y) const
   {
      PANIC("An operator cannot be used as a divisor");
      return element_type(complex());
   }

   element_type operator()(MPOperator const& x, complex const& y) const
   {
      return element_type((1.0 / y) * x);
   }
};

struct binary_power : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return std::pow(x,y);
   }

   template <typename T>
   element_type operator()(T& x, MPOperator const& y) const
   {
      PANIC("Cannot evaluate exponent of an MPOperator");
      return element_type(complex(0));
   }

   element_type operator()(MPOperator const& x, complex const& y) const
   {
      int i = int(y.real()+0.5);
      CHECK(std::norm(complex(double(i)) - y) < 1E-7)("cannot take fractional power of an operator");
      CHECK(i >= 0)("cannot take negative power of an operator");
      return i == 0 ? element_type(complex(1.0,0.0)) : pow(x, i);
   }
};

struct binary_commutator : boost::static_visitor<element_type>
{
   // commutator is zero unless we have two operators; make a template and overload for [operator,operator]
   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return element_type(complex(0.0));
   }

   // [operator,operator]
   element_type operator()(MPOperator const& x, MPOperator const& y) const
   {
      element_type Result = element_type(x*y - y*x);
      return Result;
   }
};

struct binary_dot_product : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(complex const& x, MPOperator const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(MPOperator const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(MPOperator const& x, MPOperator const& y) const
   {
      return dot(x,y);
   }
};

// the standard binary functions are put in a table

struct binary_funcs : symbols<binary_func_type>
{
   binary_funcs()
   {
      add
         ("dot", make_apply_binary_math(binary_dot_product()))
         ;
   }
};




//
// global constants
//

constants constants_p;
unary_funcs unary_funcs_p;
binary_funcs binary_funcs_p;

////////////////////////////////////////////////////////////////////////////
//
//  Semantic actions
//
////////////////////////////////////////////////////////////////////////////
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

struct push_unary
{
   push_unary(std::stack<unary_func_type>& s_) : s(s_) {}

   void operator()(unary_func_type f) const
   {
      s.push(f);
   }

   std::stack<unary_func_type>& s;
};

struct eval_unary
{
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

struct push_binary
{
   push_binary(std::stack<binary_func_type>& s_) : s(s_) {}

   void operator()(binary_func_type f) const
   {
      s.push(f);
   }

   std::stack<binary_func_type>& s;
};

struct eval_binary
{
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

template <typename BinaryFunc>
struct invoke_binary
{
   invoke_binary(std::stack<element_type >& eval_) : eval(eval_) {}
   invoke_binary(std::stack<element_type >& eval_, OperatorList const& OpList) 
      : eval(eval_), f(OpList) {}

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

struct negate_element : boost::static_visitor<element_type>
{
   template <typename T>
   element_type operator()(T const& x) const
   {
      return element_type(-x);
   }
};

struct do_negate
{
    do_negate(std::stack<element_type >& eval_)
    : eval(eval_) {}

    void operator()(char const*, char const*) const
    {
        element_type lhs = eval.top();
        eval.pop();
        eval.push(boost::apply_visitor(negate_element(), lhs));
    }

    std::stack<element_type >& eval;
};

struct push_operator
{
   push_operator(OperatorList const& OpList_, std::stack<element_type >& eval_)
      : OpList(OpList_), eval(eval_) {}

   void operator()(char const* str, char const* end) const
   {
      std::string OpName(str, end);

      CHECK(OpList.HasOperator(OpName))("Operator does not exist in lattice")(OpName);
      eval.push(element_type(OpList[OpName]));
   }

   OperatorList const& OpList;
   std::stack<element_type >& eval;
};


////////////////////////////////////////////////////////////////////////////
//
//  Our calculator grammar
//
////////////////////////////////////////////////////////////////////////////

struct calculator : public grammar<calculator>
{
    calculator(std::stack<element_type >& eval_, 
               std::stack<unary_func_type>& func_stack_,
               std::stack<binary_func_type>& bin_func_stack_,
               OperatorList const& OpList_)
       : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_), OpList(OpList_) {}

    template <typename ScannerT>
    struct definition
    {
        definition(calculator const& self)
        {
            real = ureal_p[push_real(self.eval)];

            imag = lexeme_d[ureal_p >> chset<>("iIjJ")][push_imag(self.eval)];

            identifier = lexeme_d[alpha_p >> *(alnum_p | '_')];

            bracket_expr = '(' >> *((anychar_p - chset<>("()")) | bracket_expr) >> ')';

            lattice_operator = (identifier >> !bracket_expr)[push_operator(self.OpList, self.eval)];

            unary_function = 
               eps_p(unary_funcs_p >> '(') 
               >>  unary_funcs_p[push_unary(self.func_stack)]
               >>  ('(' >> expression >> ')')[eval_unary(self.func_stack, self.eval)];

            binary_function = 
               eps_p(binary_funcs_p >> '(') 
               >>  binary_funcs_p[push_binary(self.bin_func_stack)]
               >>  ('(' >> expression >> ','  >> expression >> ')')
	       [eval_binary(self.bin_func_stack, self.eval)];

            commutator_bracket = 
               ('[' >> expression >> ',' >> expression >> ']')[invoke_binary<binary_commutator>(self.eval)];

            factor =
               imag
               |   real
               |   unary_function
	       |   binary_function
               |   keyword_d[constants_p[push_real(self.eval)]]
               |   commutator_bracket
               |   '(' >> expression >> ')'
               |   ('-' >> factor)[do_negate(self.eval)]
               |   ('+' >> factor)
               |   lattice_operator
               ;

            // power operator, next precedence, operates to the right
            pow_term =
               factor
                  >> *(  ('^' >> pow_term)[invoke_binary<binary_power>(self.eval)]
                      )
                      ;

            term =
                pow_term
                >> *(   ('*' >> pow_term)[invoke_binary<binary_multiplication>(self.eval)]
                    |   ('/' >> pow_term)[invoke_binary<binary_division>(self.eval)]
                    )
                    ;

            expression =
                term
                >> *(  ('+' >> term)[invoke_binary<binary_addition>(self.eval, self.OpList)]
                    |   ('-' >> term)[invoke_binary<binary_subtraction>(self.eval, self.OpList)]
                    )
                    ;
        }

        rule<ScannerT> expression, term, factor, real, imag, operator_literal, unary_function,
	   binary_function,
           bracket_expr, lattice_operator, identifier, pow_term, commutator_bracket;
        rule<ScannerT> const&
        start() const { return expression; }
    };

   mutable std::stack<element_type>& eval;
   mutable std::stack<unary_func_type>& func_stack;
   mutable std::stack<binary_func_type>& bin_func_stack;
   OperatorList const& OpList;
};

std::pair<OperatorList, MPOperator>
ParseLatticeAndOperator(std::string const& str)
{
   std::string::const_iterator Delim = std::find(str.begin(), str.end(), ':');
   if (Delim == str.end())
   {
      PANIC("fatal: expression of the form \"lattice:expression\" expected.")(str);
   }

   std::string Lattice = std::string(str.begin(), Delim);
   boost::trim(Lattice);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(Lattice);

   ++Delim;
   std::string Expr(Delim, str.end());

   std::stack<element_type> eval;
   std::stack<unary_func_type> func_stack;
   std::stack<binary_func_type> bin_func_stack;
   
   calculator calc(eval, func_stack, bin_func_stack, *System); //  Our parser
   parse_info<> info = parse(Expr.c_str(), calc, space_p);
   if (!info.full)
   {
      PANIC("Operator parser failed, stopped at")(info.stop);
   }

   CHECK(func_stack.empty());
   CHECK(!eval.empty());
   element_type Result = eval.top();
   eval.pop();
   CHECK(eval.empty());

   MPOperator* Op = boost::get<MPOperator>(&Result);
   if (Op)
   {
      return std::make_pair(*System, *Op);
   }
   // else
   complex x = boost::get<complex>(Result);

   return std::make_pair(*System, x * (*System)["I"]);
}
