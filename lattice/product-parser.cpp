// -*- C++ -*- $Id$

#include "product-parser.h"
#include "parser/parser.h"
#include "unitcell-parser.h"
#include <boost/algorithm/string.hpp>


using namespace Parser;

typedef boost::variant<complex, ProductMPO> element_type;

// specializations of some templates
namespace Parser
{

template <>
struct ElementExp<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return std::exp(c);
   }

   element_type operator()(ProductMPO const&) const
   {
      PANIC("exp() is not implemented for ProductMPO's!");
      return ProductMPO();
   }
};

template <>
struct negate_element<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex c) const
   {
      return -c;
   }

   template <typename T>
   element_type operator()(T const& x) const
   {
      PANIC("Unary negation is not defined for ProductOperator");
      return element_type(0.0);
   }
};


template <>
struct binary_power<element_type> : boost::static_visitor<element_type>
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
   element_type operator()(complex const&, ProductMPO const&) const
   {
      PANIC("exp() is not implemented for ProductMPO's!");
      return ProductMPO();
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Cannot evaluate exponent with an operator")(typeid(y).name());
      return element_type(complex(0));
   }
};

template <>
struct binary_addition<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   { 
      return element_type(x+y);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Addition is not defined for ProductOperator");
      return complex(0.0);
   }
};

template <>
struct binary_subtraction<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x-y);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Subtraction is not defined for ProductOperator");
      return complex(0.0);
   }
};

template <>
struct binary_multiplication<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      return x*y;
   }
};

template <>
struct binary_division<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x/y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      PANIC("Cannot divide a scalar by a ProductOperator");
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      PANIC("Cannot divide a ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Cannot divide ProductOperators");
      return complex(0.0);
   }
};

template <>
struct binary_dot_product<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      return dot(x,y);
   }
};

template <>
struct binary_inner_product<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      return inner(x,y);
   }
};

template <>
struct binary_outer_product<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   template <typename T>
   element_type operator()(complex const& x, T const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T>
   element_type operator()(T const& x, complex const& y) const
   {
      PANIC("Cannot multiply ProductOperator by a scalar");
      return complex(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      return x*y;
   }
};

template <>
struct binary_cross_product<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Cross product is not defined for ProductOperator");
      return element_type(0.0);
   }
};

template <>
struct binary_commutator<element_type> : boost::static_visitor<element_type>
{
   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(0.0);
   }

   template <typename T, typename U>
   element_type operator()(T const& x, U const& y) const
   {
      PANIC("Commutator is not defined for ProductOperator");
      return element_type(0.0);
   }
};

} // namespace Parser

namespace ILP // to avoid confusion of duplicate names
{

struct push_operator
{
   push_operator(InfiniteLattice const& Lattice_,
		 std::stack<std::string>& IdentStack_, 
		 std::stack<element_type>& eval_)
      : Lattice(Lattice_), IdentStack(IdentStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      CHECK(Lattice.product_operator_exists(OpName))("Operator does not exist in the lattice")(OpName);

      eval.push(element_type(Lattice.Product(OpName)));
   }

   InfiniteLattice const& Lattice;
   std::stack<std::string>& IdentStack;
   std::stack<element_type >& eval;
};

struct push_operator_param
{
   push_operator_param(InfiniteLattice const& Lattice_,
		       std::stack<std::string>& IdentStack_, 
		       std::stack<std::stack<element_type> >& ParamStack_,
		       std::stack<element_type>& eval_)
      : Lattice(Lattice_), IdentStack(IdentStack_), 
	ParamStack(ParamStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      // Convert the parameter stack into a vector of numbers
      int NumParam = ParamStack.top().size(); // number of parameters
      std::vector<std::complex<double> > Params(NumParam);
      while (!ParamStack.top().empty())
      {
	 Params[--NumParam] = boost::get<std::complex<double> >(ParamStack.top().top());
	 ParamStack.top().pop();
      }
      ParamStack.pop();

      eval.push(element_type(Lattice.ProductOperatorFunction(OpName, Params)));
   }

   InfiniteLattice const& Lattice;
   std::stack<std::string>& IdentStack;
   std::stack<std::stack<element_type> >& ParamStack;
   std::stack<element_type >& eval;
};

// We read parameters initially into the element stack, so 
// this function moves an item from the element stack to the parameter stack
template <typename T>
struct move_to_parameter_pack;

template <>
struct move_to_parameter_pack<element_type>
{
   move_to_parameter_pack(std::stack<std::stack<element_type> >& ParamStack_,
			  std::stack<element_type>& eval_)
      : ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      ParamStack.top().push(eval.top());
      eval.pop();
   }

   std::stack<std::stack<element_type> >& ParamStack;
   std::stack<element_type >& eval;
};

template <typename T>
struct push_parameter_pack;

template <>
struct push_parameter_pack<element_type>
{
   push_parameter_pack(std::stack<std::stack<element_type> >& ParamStack_)
      : ParamStack(ParamStack_) {}

   void operator()(char const*, char const*) const
   {
      ParamStack.push(std::stack<element_type>());
   }

   void operator()(char) const
   {
      ParamStack.push(std::stack<element_type>());
   }

   std::stack<std::stack<element_type> >& ParamStack;
};

struct push_prod_unit
{
   push_prod_unit(InfiniteLattice const& Lattice_,
		 std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));

      eval.push(prod_unit_left_to_right(Op.MPO(), Lattice.GetUnitCell().size()));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

struct push_prod_unit_r
{
   push_prod_unit_r(InfiniteLattice const& Lattice_,
		    std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));

      eval.push(prod_unit_right_to_left(Op.MPO(), Lattice.GetUnitCell().size()));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

struct push_string
{
   push_string(InfiniteLattice const& Lattice_, std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      FiniteMPO Op = ParseStringOperator(*Lattice.GetUnitCell().GetSiteList(), 
					 std::string(Start, End), Lattice.GetUnitCell().size());
      eval.push(prod_unit_left_to_right(Op, Lattice.GetUnitCell().size()));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

} // namespace ILP;

using namespace ILP;

struct ProductParser : public grammar<ProductParser>
{
   typedef boost::variant<complex, ProductMPO> element_type;
   typedef boost::function<element_type(element_type)> unary_func_type;
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;

   typedef std::stack<element_type>       ElemStackType;
   typedef std::stack<std::string>        IdentStackType;
   typedef std::stack<unary_func_type>    UnaryFuncStackType;
   typedef std::stack<binary_func_type>   BinaryFuncStackType;
   typedef std::stack<ElemStackType>      ParamStackType;

   static constants                  constants_p;
   static unary_funcs<element_type>  unary_funcs_p;
   //   static binary_funcs<element_type> binary_funcs_p;

   ProductParser(ElemStackType& eval_, 
		 IdentStackType& identifier_stack_,
		 ParamStackType& param_stack_,
		 UnaryFuncStackType& func_stack_,
		 BinaryFuncStackType& bin_func_stack_,
		 InfiniteLattice const& Lattice_)
      : eval(eval_), identifier_stack(identifier_stack_), 
	param_stack(param_stack_),
	func_stack(func_stack_), bin_func_stack(bin_func_stack_), Lattice(Lattice_),
	binary_funcs_p() {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(ProductParser const& self)
      {
	 real = ureal_p[push_real<element_type>(self.eval)];
	 
	 imag = lexeme_d[ureal_p >> chset<>("iIjJ")][push_imag<element_type>(self.eval)];
	 
	 identifier = lexeme_d[alpha_p >> *(alnum_p | '_')][push_identifier(self.identifier_stack)];

	 // We re-use the identifier_stack for quantum numbers, and rely on the grammar rules to
	 // avoid chaos!
	 quantumnumber = lexeme_d[*(anychar_p - chset<>("()"))]
	    [push_identifier(self.identifier_stack)];

	 parameter = expression[move_to_parameter_pack<element_type>(self.param_stack, self.eval)];

	 parameter_list = ch_p('{')[push_parameter_pack<element_type>(self.param_stack)]
	    >> parameter >> *(',' >> parameter) >> '}';

	 bracket_expr = '(' >> expression >> ')';

	 sq_bracket_expr = '[' >> expression >> ']';

	 expression_string = lexeme_d[+((anychar_p - chset<>("()"))
					| (ch_p('(') >> expression_string >> ch_p(')')))];

	 string_expression = str_p("string")
	    >> '(' 
	    >> expression_string[push_string(self.Lattice, self.eval)]
	    >> ')';

	 prod_unit_expression = str_p("prod_unit")
	    >> '(' 
	    >> expression_string[push_prod_unit(self.Lattice, self.eval)]
	    >> ')';

	 prod_unit_r_expression = str_p("prod_unit_r")
	    >> '(' 
	    >> expression_string[push_prod_unit_r(self.Lattice, self.eval)]
	    >> ')';

	 operator_expression = 
	    identifier
	    >> (parameter_list[push_operator_param(self.Lattice, 
						   self.identifier_stack, 
						   self.param_stack,
						   self.eval)]
		|   eps_p[push_operator(self.Lattice, self.identifier_stack, self.eval)]
		);
	 
	 unary_function = 
	    eps_p(unary_funcs_p >> '(') 
	    >>  unary_funcs_p[push_unary<element_type>(self.func_stack)]
	    >>  ('(' >> expression >> ')')[eval_unary<element_type>(self.func_stack, self.eval)];
	 
	 binary_function = 
	    eps_p(self.binary_funcs_p >> '(') 
	    >>  self.binary_funcs_p[push_binary<element_type>(self.bin_func_stack)]
	    >>  ('(' >> expression >> ','  >> expression >> ')')
	    [eval_binary<element_type>(self.bin_func_stack, self.eval)];
	 
	 commutator_bracket = 
	    ('[' >> expression >> ',' >> expression >> ']')[invoke_binary<element_type, 
							    binary_commutator<element_type> >(self.eval)];
	 
	 factor =
	    imag
	    |   real
	    |   unary_function
	    |   binary_function
	    |   keyword_d[constants_p[push_real<element_type>(self.eval)]]
	    |   prod_unit_expression
	    |   prod_unit_r_expression
	    |   string_expression
	    |   '(' >> expression >> ')'
	    |   ('-' >> factor)[do_negate<element_type>(self.eval)]
	    |   ('+' >> factor)
	    |   operator_expression
	    ;
	 
	 // power operator, next precedence, operates to the right
	 pow_term =
	    factor
	    >> *(  ('^' >> pow_term)[invoke_binary<element_type, binary_power<element_type> >(self.eval)]
		   )
	    ;
	 
	 term =
	    pow_term
	    >> *(   ('*' >> pow_term)[invoke_binary<element_type, 
				      binary_multiplication<element_type> >(self.eval)]
                    |   ('/' >> pow_term)[invoke_binary<element_type, 
					  binary_division<element_type> >(self.eval)]
                    )
	    ;
	 
	 expression =
	    term
	    >> *(  ('+' >> term)[invoke_binary<element_type, 
				 binary_addition<element_type> >(self.eval)]
		   |   ('-' >> term)[invoke_binary<element_type, 
				     binary_subtraction<element_type> >(self.eval)]
		   )
	    ;
      }
      
      rule<ScannerT> expression, term, factor, real, imag, operator_literal, unary_function,
	 binary_function, bracket_expr, quantumnumber, sq_bracket_expr, 
	 operator_expression, operator_bracket_sq, operator_sq_bracket, operator_bracket, operator_sq,
	 parameter, parameter_list, expression_string, prod_unit_expression, prod_unit_r_expression,
	 string_expression,
	 identifier, pow_term, commutator_bracket;
      rule<ScannerT> const&
      start() const { return expression; }
   };
   
   std::stack<element_type>& eval;
   IdentStackType& identifier_stack;
   ParamStackType& param_stack;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   InfiniteLattice const& Lattice;

   binary_funcs<element_type> binary_funcs_p;
};


// global variables (static members of ProductParser)
constants ProductParser::constants_p;
unary_funcs<ProductParser::element_type> ProductParser::unary_funcs_p;
//binary_funcs<ProductParser::element_type> ProductParser::binary_funcs_p;

ProductMPO
ParseProductOperator(InfiniteLattice const& Lattice, std::string const& Str)
{
   typedef ProductParser::element_type element_type;

   ProductParser::ElemStackType ElemStack;
   ProductParser::IdentStackType IdentStack;
   ProductParser::ParamStackType ParamStack;
   ProductParser::UnaryFuncStackType UnaryFuncStack;
   ProductParser::BinaryFuncStackType BinaryFuncStack;

   ProductParser Parser(ElemStack, IdentStack, ParamStack, 
			UnaryFuncStack, BinaryFuncStack, Lattice);

   parse_info<> info = parse(Str.c_str(), Parser, space_p);
   if (!info.full)
   {
      PANIC("Operator parser failed, stopped at")(info.stop);
   }

   CHECK(UnaryFuncStack.empty());
   CHECK(BinaryFuncStack.empty());
   CHECK(IdentStack.empty());
   CHECK(ParamStack.empty());
   CHECK(!ElemStack.empty());
   element_type Result = ElemStack.top();
   ElemStack.pop();
   CHECK(ElemStack.empty());

   ProductMPO* Op = boost::get<ProductMPO>(&Result);
   if (Op)
   {
      return *Op;
   }
   // we cannot handle a pure number
   PANIC("Cannot convert a number into a ProductMPO");
   return ProductMPO();
}

std::pair<ProductMPO, InfiniteLattice>
ParseProductOperatorAndLattice(std::string const& Str)
{
   std::string::const_iterator Delim = std::find(Str.begin(), Str.end(), ':');
   if (Delim == Str.end())
   {
      PANIC("fatal: expression of the form \"lattice:expression\" expected.")(Str);
   }

   std::string LatticeFile = std::string(Str.begin(), Delim);
   boost::trim(LatticeFile);
   pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

   ++Delim;
   std::string Expr(Delim, Str.end());

   ProductMPO Op = ParseProductOperator(*Lattice, Expr);
   return std::make_pair(Op, *Lattice);
}
