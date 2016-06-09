// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/triangular-parser.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "triangular-parser.h"
#include "parser/parser.h"
#include "unitcell-parser.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>

using namespace Parser;

typedef boost::variant<complex, TriangularMPO> element_type;

namespace ILP // to avoid confusion of duplicate names
{

// utility to pop an integer off the element stack
int pop_int(std::stack<element_type>& eval)
{
   complex x = boost::get<complex>(eval.top());
   eval.pop();
   int j = boost::math::iround(x.real());
   CHECK(norm_frob(x - double(j)) < 1E-7)("index must be an integer")(x);
   return j;
}

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

      CHECK(Lattice.triangular_operator_exists(OpName))("Operator does not exist in the lattice")(OpName);

      eval.push(element_type(Lattice.Triangular(OpName)));
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

      eval.push(element_type(Lattice.TriangularOperatorFunction(OpName, Params)));
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

// convert the top of the expression stack from a number of cells
// to a number of sites, by multiplying by the unit cell size
struct scale_cells_to_sites
{
   scale_cells_to_sites(InfiniteLattice const& Lattice_,
			std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      int Cells = pop_int(eval);
      eval.push(std::complex<double>(double(Cells*Lattice.GetUnitCell().size())));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

// push the number of sites in the lattice onto the top of the expression stack
struct push_number_of_sites
{
   push_number_of_sites(InfiniteLattice const& Lattice_,
			std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      eval.push(std::complex<double>(double(Lattice.GetUnitCell().size())));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

struct push_sum_unit
{
   push_sum_unit(InfiniteLattice const& Lattice_,
		 std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End))(Sites);
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(sum_unit(Op, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

struct push_sum_k
{
   push_sum_k(InfiniteLattice const& Lattice_,
	      std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      std::complex<double> k = boost::get<std::complex<double> >(eval.top());
      eval.pop();
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(sum_k(k, Op, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

struct push_sum_kink
{
   push_sum_kink(InfiniteLattice const& Lattice_,
		 std::stack<element_type>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();

      // Find the comma separating the two operators
      int nBracket = 0;
      char const* Comma = Start;
      while (Comma != End && (*Comma != ',' || nBracket > 0))
      {
	 if (*Comma == '(')
	    ++nBracket;
	 else if (*Comma == ')')
	    --nBracket;
	 ++Comma;
      }
      if (Comma == End)
      {
	 PANIC("Failed to parse two parameters in sum_kink");
      }

      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Kink = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma));
      ++Comma; // skip over the comma
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Comma, End));
      Op.ExtendToCoverUnitCell(Sites);

      eval.push(sum_kink(Kink, Op, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<element_type >& eval;
};

} // namespace ILP;

using namespace ILP;

struct InfiniteLatticeParser : public grammar<InfiniteLatticeParser>
{
   typedef boost::variant<complex, TriangularMPO> element_type;
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

   InfiniteLatticeParser(ElemStackType& eval_, 
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
      definition(InfiniteLatticeParser const& self)
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

	 prod_expression = (str_p("prod") >> '(' >> expression >> ',' >> expression >> ',' >> quantumnumber >> ')')
	    [push_prod<element_type>(self.identifier_stack, self.eval)];
 
	 bracket_expr = '(' >> expression >> ')';

	 sq_bracket_expr = '[' >> expression >> ']';

	 expression_string = lexeme_d[+((anychar_p - chset<>("()"))
					| (ch_p('(') >> expression_string >> ch_p(')')))];

	 num_cells = (eps_p((str_p("cells") | str_p("sites")) >> '=')
		      >> ((str_p("cells") >> '=' >> expression >> ',')
			  [scale_cells_to_sites(self.Lattice, self.eval)]
			  | (str_p("sites") >> '=' >> expression >> ',')))
	    | eps_p[push_number_of_sites(self.Lattice, self.eval)];
      
	 sum_unit_expression = str_p("sum_unit")
	    >> '('
	    >> num_cells
	    >> expression_string[push_sum_unit(self.Lattice, self.eval)]
	    >> ')';

	 sum_kink_expression = str_p("sum_kink")
	    >> '(' 
	    >> num_cells
	    >> expression_string[push_sum_kink(self.Lattice, self.eval)]
	    >> ')';
	 
	 sum_k_expression = str_p("sum_k")
	    >> '(' 
	    >> num_cells
	    >> expression >> ','
	    >> expression_string[push_sum_k(self.Lattice, self.eval)]
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
	    |   prod_expression
	    |   sum_unit_expression
	    |   sum_kink_expression
	    |   sum_k_expression
	    |   commutator_bracket
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
	 binary_function, bracket_expr, quantumnumber, prod_expression, sq_bracket_expr, 
	 operator_expression, operator_bracket_sq, operator_sq_bracket, operator_bracket, operator_sq,
	 parameter, parameter_list, expression_string, sum_unit_expression, sum_kink_expression, sum_k_expression,
	 identifier, pow_term, commutator_bracket, num_cells;
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


// global variables (static members of InfiniteLatticeParser)
constants InfiniteLatticeParser::constants_p;
unary_funcs<InfiniteLatticeParser::element_type> InfiniteLatticeParser::unary_funcs_p;
//binary_funcs<InfiniteLatticeParser::element_type> InfiniteLatticeParser::binary_funcs_p;

TriangularMPO
ParseTriangularOperator(InfiniteLattice const& Lattice, std::string const& Str)
{
   typedef InfiniteLatticeParser::element_type element_type;

   InfiniteLatticeParser::ElemStackType ElemStack;
   InfiniteLatticeParser::IdentStackType IdentStack;
   InfiniteLatticeParser::ParamStackType ParamStack;
   InfiniteLatticeParser::UnaryFuncStackType UnaryFuncStack;
   InfiniteLatticeParser::BinaryFuncStackType BinaryFuncStack;

   InfiniteLatticeParser Parser(ElemStack, IdentStack, ParamStack, 
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

   TriangularMPO* Op = boost::get<TriangularMPO>(&Result);
   if (Op)
   {
      return *Op;
   }
   // else, we also handle the case where the operator is a number
   complex x = boost::get<complex>(Result);
   return x*Lattice.Triangular("I");
}

std::pair<TriangularMPO, InfiniteLattice>
ParseTriangularOperatorAndLattice(std::string const& Str)
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

   TriangularMPO Op = ParseTriangularOperator(*Lattice, Expr);
   return std::make_pair(Op, *Lattice);
}
