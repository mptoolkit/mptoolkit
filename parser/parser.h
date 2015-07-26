// -*- C++ -*- $Id$
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

#include <boost/spirit/include/classic_core.hpp>
#include <boost/spirit/include/classic_chset.hpp>
#include <boost/spirit/include/classic_symbols.hpp>
#include <boost/spirit/include/classic_refactoring.hpp>
#include <boost/spirit/include/classic_distinct.hpp>

//#include <boost/spirit/include/qi_core.hpp>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/variant.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/get.hpp>
#include <boost/math/complex.hpp>
#include <iostream>
#include <stack>
#include <functional>
#include <string>
#include <complex>

namespace Parser
{

inline
std::string Spaces(int Size)
{
   return std::string(Size, ' ');
}

class ParserError : public std::exception
{
   public:
      explicit ParserError(std::string const& Why);

      ParserError(ParserError const& Prev, std::string const& Why);

      ParserError(std::exception const& Prev, std::string const& Why);

      ~ParserError() throw() { }

      virtual const char* what() const throw() { return Msg.c_str(); }

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
      ParserError(std::list<std::string> const& CallStack_, char const* Position);
      ParserError(std::list<std::string> const& CallStack_, char const* Position, char const* End_);
      ParserError(std::list<std::string> const& CallStack_, 
		  std::string const& Why, char const* Position, char const* End_,
		  char const* beg, char const* end);

      void AssembleMessage();

      std::list<std::string> CallStack;
      std::string Msg;

      char const* Pos;
      char const* End;
};

inline
ParserError::ParserError(std::string const& Why)
   : CallStack(1, Why), Pos(NULL), End(NULL)
{
   this->AssembleMessage();
}

inline
ParserError::ParserError(ParserError const& Prev, std::string const& Why)
   : CallStack(Prev.CallStack), Pos(Prev.Pos), End(Prev.End)
{
   CallStack.push_back(Why);
   this->AssembleMessage();
}

inline
ParserError::ParserError(std::exception const& Prev, std::string const& Why)
   : CallStack(1, Prev.what()), Pos(NULL), End(NULL)
{
   CallStack.push_back(Why);
   this->AssembleMessage();
}


inline
ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position)
   : CallStack(CallStack_), Pos(Position), End(NULL)
{
   this->AssembleMessage();
}

inline
ParserError::ParserError(std::list<std::string> const& CallStack_, char const* Position,
			 char const* End_)
   : CallStack(CallStack_), Pos(Position), End(End_)
{
   this->AssembleMessage();
}

inline
ParserError::ParserError(std::list<std::string> const& CallStack_, 
			 std::string const& Why, char const* Position, char const* End_,
			 char const* beg, char const* end)
   : CallStack(CallStack_), Pos(NULL), End(NULL)
{
   std::string Next = Why + "\n" + std::string(beg, end);
   if (Position)
   {
      int Count = End_ ? std::distance(Position, End_) : 1;
      Next = Next + "\n" + std::string(std::distance(beg, Position), ' ')
	 + std::string(Count, '^');
   }
   CallStack.push_back(Next);
   this->AssembleMessage();
}

inline
ParserError 
ParserError::AtPosition(std::string const& Why, char const* Position)
{
   return ParserError(std::list<std::string>(1, Why), Position);
}

inline
ParserError 
ParserError::AtPosition(ParserError const& Prev, char const* Position)
{
   return ParserError(Prev.CallStack, Position);
}

inline
ParserError 
ParserError::AtPosition(std::exception const& Prev, char const* Position)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Position);
}

inline
ParserError
ParserError::AtRange(std::string const& Why, char const* Start, char const* End)
{
   return ParserError(std::list<std::string>(1, Why), Start, End);
}

inline
ParserError
ParserError::AtRange(ParserError const& Prev, char const* Start, char const* End)
{
   return ParserError(Prev.CallStack, Start, End);
}

inline
ParserError
ParserError::AtRange(std::exception const& Prev, char const* Start, char const* End)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Start, End);
}

inline
ParserError
ParserError::Finalize(ParserError const& Prev, std::string const& Why, char const* beg, char const* end)
{
   return ParserError(Prev.CallStack, Why, Prev.Pos, Prev.End, beg, end);
}

inline
ParserError
ParserError::Finalize(std::exception const& Prev, std::string const& Why, char const* beg, char const* end)
{
   return ParserError(std::list<std::string>(1, Prev.what()), Why, NULL, NULL, beg, end);
}

inline
void
ParserError::AssembleMessage()
{
   Msg = "Parser error:\n";
   for (std::list<std::string>::const_iterator I = CallStack.begin(); I != CallStack.end(); ++I)
   {
      if (I != CallStack.begin())
	 Msg = Msg + '\n';
      Msg = Msg + (*I);
   }
}

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
	    s = s + "\n" + Spaces(std::distance(beg, I)) + "^";
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

using namespace boost::spirit::classic;

distinct_parser<> const keyword_p("0-9a-zA-Z_");

distinct_directive<> const keyword_d("0-9a-zA-Z_");

typedef std::complex<double> complex;

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
   return boost::apply_visitor(NameOf(), x);
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
   if (!boost::get<complex>(&eval.top()))
      throw ParserError("expected an integer, got a " + name_of(eval.top()));
   complex x = boost::get<complex>(eval.top());
   eval.pop();
   int j = boost::math::iround(x.real());
   if (LinearAlgebra::norm_frob(x - double(j)) > 1E-7)
       throw ParserError("expected an integer, got a real/complex number: " + format_complex(x));
   return j;
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
   return boost::math::asin(x);
}

inline
complex cacos(complex x)
{
   return boost::math::acos(x);
}

inline
complex catan(complex x)
{
   return boost::math::atan(x);
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
struct unary_funcs : symbols<boost::function<element_type(element_type)> >
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

   void operator()(char const* a, char const* b) const
   {
      try
      {
	 element_type rhs = eval.top();
	 eval.pop();
	 element_type lhs = eval.top();
	 eval.pop();
	 eval.push(boost::apply_visitor(f, lhs, rhs));
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
	 eval.push(boost::apply_visitor(negate_element<element_type>(), lhs));
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
	 
	 eval.push(boost::apply_visitor(ternary_product<element_type>(q), op1, op2));
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
      ParamStack.top().push_back(Function::Parameter(boost::get<complex>(eval.top())));
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
						     boost::get<complex>(eval.top())));
      IdentifierStack.pop();
      eval.pop();
   }

   std::stack<element_type>& eval;
   std::stack<std::string>& IdentifierStack;
   std::stack<Function::ParameterList>& ParamStack;
};

} //namespace Parser

#endif
