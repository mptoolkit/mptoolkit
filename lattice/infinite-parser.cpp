// -*- C++ -*- $Id: triangular-parser.cpp 1482 2015-05-13 05:47:05Z ianmcc $

#include "infinite-parser.h"
#include "mpo/infinite_mpo_actions.h"
#include "parser/parser.h"
#include "unitcell-parser.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>

//
// TODO: we can allow invoking complex-valued unit cell functions and operator functions,
//

namespace ILP // to avoid confusion of duplicate names
{

typedef InfiniteMPOElement ElementType;

struct push_operator
{
   push_operator(InfiniteLattice const& Lattice_,
		 std::stack<std::string>& IdentStack_, 
		 std::stack<ElementType>& eval_)
      : Lattice(Lattice_), IdentStack(IdentStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();
      CHECK(Lattice.operator_exists(OpName))("Operator does not exist in the lattice")(OpName);
      eval.push(Lattice[OpName].op());
   }

   InfiniteLattice const& Lattice;
   std::stack<std::string>& IdentStack;
   std::stack<ElementType>& eval;
};

struct eval_function
{
   eval_function(InfiniteLattice const& Lattice_,
		 std::stack<std::string>& FunctionStack_, 
		 std::stack<Function::ParameterList>& ParamStack_,
		 std::stack<ElementType>& eval_)
      : Lattice(Lattice_), FunctionStack(FunctionStack_), 
	ParamStack(ParamStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      eval.push(Lattice.eval_function(FunctionStack.top(), ParamStack.top()).op());
      FunctionStack.pop();
      ParamStack.pop();
   }

   InfiniteLattice const& Lattice;
   std::stack<std::string>& FunctionStack;
   std::stack<Function::ParameterList>& ParamStack;
   std::stack<ElementType>& eval;
};

// convert the top of the expression stack from a number of cells
// to a number of sites, by multiplying by the unit cell size
struct scale_cells_to_sites
{
   scale_cells_to_sites(InfiniteLattice const& Lattice_,
			std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* a, char const*) const
   {
      try
      {
	 int Cells = pop_int(eval);
	 eval.push(std::complex<double>(double(Cells*Lattice.GetUnitCell().size())));
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

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

// push the number of sites in the lattice onto the top of the expression stack
struct push_number_of_sites
{
   push_number_of_sites(InfiniteLattice const& Lattice_,
			std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      eval.push(std::complex<double>(double(Lattice.GetUnitCell().size())));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_prod_unit
{
   push_prod_unit(InfiniteLattice const& Lattice_,
		 std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(prod_unit_left_to_right(Op.MPO(), Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_prod_unit_r
{
   push_prod_unit_r(InfiniteLattice const& Lattice_,
		    std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(prod_unit_right_to_left(Op.MPO(), Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_string
{
   push_string(InfiniteLattice const& Lattice_, std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      FiniteMPO Op = ParseStringOperator(*Lattice.GetUnitCell().GetSiteList(), 
					 std::string(Start, End), Lattice.GetUnitCell().size());
      eval.push(prod_unit_left_to_right(Op, Lattice.GetUnitCell().size()));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_sum_unit
{
   push_sum_unit(InfiniteLattice const& Lattice_,
		 std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      try
      {
	 int Sites = pop_int(eval);
	 int Cells = Sites / Lattice.GetUnitCell().size();
	 DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End))(Sites);
	 UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));
	 //Op.ExtendToCoverUnitCell(Sites);
	 eval.push(sum_unit(Op, Sites));
      }
      catch (ParserError const& p)
      {
	 throw ParserError::AtPosition(p, Start);
      }
      catch (std::exception const& p)
      {
	 throw ParserError::AtPosition(p, Start);
      }
      catch (...)
      {
	 throw;
      }
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_sum_k
{
   push_sum_k(InfiniteLattice const& Lattice_,
	      std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      std::complex<double> k = boost::get<std::complex<double> >(eval.top());
      eval.pop();
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End));
      //Op.ExtendToCoverUnitCell(Sites);
      eval.push(sum_k(k, Op, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_sum_kink
{
   push_sum_kink(InfiniteLattice const& Lattice_,
		 std::stack<ElementType>& eval_)
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
      //Op.ExtendToCoverUnitCell(Sites);
      eval.push(sum_kink(Kink, Op, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

struct push_sum_string_inner
{
   push_sum_string_inner(InfiniteLattice const& Lattice_,
		 std::stack<ElementType>& eval_)
      : Lattice(Lattice_), eval(eval_) {}
   
   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      int Cells = Sites / Lattice.GetUnitCell().size();

      // here we expect 3 operators, separated by commas
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
	 PANIC("Failed to parse three parameters in sum_string_inner");
      }

      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op1 = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma));
      // find the next comma
      ++Comma;
      Start = Comma;
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
	 PANIC("Failed to parse three parameters in sum_string_inner");
      }
      UnitCellMPO String = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma));
      ++Comma; // skip over the comma
      UnitCellMPO Op2 = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Comma, End));
      eval.push(sum_string_inner(Op1, String, Op2, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
};

} // namespace ILP;

using namespace ILP;

struct InfiniteLatticeParser : public grammar<InfiniteLatticeParser>
{
   typedef InfiniteMPOElement ElementType;
   typedef boost::function<ElementType(ElementType)> unary_func_type;
   typedef boost::function<ElementType(ElementType, ElementType)> binary_func_type;

   typedef std::stack<ElementType>             ElemStackType;
   typedef std::stack<unary_func_type>         UnaryFuncStackType;
   typedef std::stack<binary_func_type>        BinaryFuncStackType;
   typedef std::stack<std::string>             FunctionStackType;
   typedef std::stack<std::string>             IdentifierStackType;
   typedef std::stack<Function::ParameterList> ParameterStackType;
   typedef symbols<complex>                    ArgumentType;

   static constants                  constants_p;
   static unary_funcs<ElementType>  unary_funcs_p;
   static binary_funcs<ElementType> binary_funcs_p;

   InfiniteLatticeParser(ElemStackType& eval_, 
		  UnaryFuncStackType& func_stack_,
		  BinaryFuncStackType& bin_func_stack_,
		  IdentifierStackType& IdentifierStack_,
		  FunctionStackType& Functions_,
		  ParameterStackType& Parameters_,
		  ArgumentType& Arguments_,
		  InfiniteLattice const& Lattice_)
      : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_),
      IdentifierStack(IdentifierStack_), FunctionStack(Functions_),
      ParameterStack(Parameters_), Arguments(Arguments_),
	Lattice(Lattice_)
   {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(InfiniteLatticeParser const& self)
      {
	 real = ureal_p[push_value<ElementType>(self.eval)];
	 
	 imag = lexeme_d[ureal_p >> chset<>("iIjJ")][push_imag<ElementType>(self.eval)];
	 
	 identifier = lexeme_d[alpha_p >> *(alnum_p | '_')];
	 //[push_identifier(self.IdentifierStack)];

	 // We re-use the IdentifierStack for quantum numbers, and rely on the grammar rules to
	 // avoid chaos!
	 quantumnumber = lexeme_d[*(anychar_p - chset<>("()"))]
	    [push_identifier(self.IdentifierStack)];

	 named_parameter = eps_p(identifier >> '=')
	    >> identifier[push_identifier(self.IdentifierStack)]
	    >> '='
	    >> expression[push_named_parameter<ElementType>(self.eval, 
							     self.IdentifierStack, 
							     self.ParameterStack)];

	 parameter = expression[push_parameter<ElementType>(self.eval, self.ParameterStack)];

	 // parameter_list is a comma-separated list of parameters, may be empty
	 // at least one parameter
	 parameter_list = '{' >> (!((named_parameter | parameter) % ',')) >> '}';

	 prod_expression = (str_p("prod") >> '(' >> expression >> ',' >> expression >> ',' >> quantumnumber >> ')')
	    [push_prod<ElementType>(self.IdentifierStack, self.eval)];
 
	 bracket_expr = '(' >> expression >> ')';

	 sq_bracket_expr = '[' >> expression >> ']';

	 expression_string = lexeme_d[+((anychar_p - chset<>("()"))
					| (ch_p('(') >> expression_string >> ch_p(')')))];

	 num_cells = (eps_p((str_p("cells") | str_p("sites")) >> '=')
		      >> ((str_p("cells") >> '=' >> expression >> ',')
			  [scale_cells_to_sites(self.Lattice, self.eval)]
			  | (str_p("sites") >> '=' >> expression >> ',')))
	    | eps_p[push_number_of_sites(self.Lattice, self.eval)];
      
	 // ProductMPO expressions

	 string_expression = str_p("string")
	    >> '(' 
	    >> expression_string[push_string(self.Lattice, self.eval)]
	    >> ')';

	 prod_unit_expression = str_p("prod_unit")
	    >> '(' 
	    >> num_cells
	    >> expression_string[push_prod_unit(self.Lattice, self.eval)]
	    >> ')';

	 prod_unit_r_expression = str_p("prod_unit_r")
	    >> '(' 
	    >> num_cells
	    >> expression_string[push_prod_unit_r(self.Lattice, self.eval)]
	    >> ')';

	 // TriangularMPO expressions

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

	 sum_string_inner_expression = str_p("sum_string_inner")
	    >> '('
	    >> num_cells
	    >> expression_string[push_sum_string_inner(self.Lattice, self.eval)]
	    >> ')';

#if 0
	 sum_string_dot_expression = str_p("sum_string_dot")
	    >> '('
	    >> num_cells
	    >> expression_string[push_sum_string_dot(self.Lattice, self.eval)]
	    >> ')';
#endif

	 function_expression = eps_p(identifier >> '{')
	    >> identifier[push_function(self.FunctionStack, self.ParameterStack)]
	    >> parameter_list[eval_function(self.Lattice, 
					     self.FunctionStack, 
					     self.ParameterStack,
					     self.eval)];
	 
	 operator_expression = 
		identifier[push_identifier(self.IdentifierStack)]
		[push_operator(self.Lattice, self.IdentifierStack, self.eval)];
	 
	 unary_function = 
	    eps_p(unary_funcs_p >> '(') 
	    >>  unary_funcs_p[push_unary<ElementType>(self.func_stack)]
	    >>  ('(' >> expression >> ')')[eval_unary<ElementType>(self.func_stack, self.eval)];
	 
	 binary_function = 
	    eps_p(self.binary_funcs_p >> '(') 
	    >>  self.binary_funcs_p[push_binary<ElementType>(self.bin_func_stack)]
	    >>  ('(' >> expression >> ','  >> expression >> ')')
	    [eval_binary<ElementType>(self.bin_func_stack, self.eval)];
	 
	 commutator_bracket = 
	    ('[' >> expression >> ',' >> expression >> ']')[invoke_binary<ElementType, 
							    binary_commutator<ElementType> >(self.eval)];
	 
	 factor =
	    imag
	    |   real
	    |   unary_function
	    |   binary_function
	    |   keyword_d[constants_p[push_value<ElementType>(self.eval)]]
	    |   keyword_d[self.Arguments[push_value<ElementType>(self.eval)]]
	    |   prod_expression
	    |   prod_unit_expression
	    |   prod_unit_r_expression
	    |   string_expression
	    |   sum_unit_expression
	    |   sum_kink_expression
	    |   sum_k_expression
	    |   sum_string_inner_expression
	    |   sum_string_dot_expression
	    |   commutator_bracket
	    |   '(' >> expression >> ')'
	    |   ('-' >> factor)[do_negate<ElementType>(self.eval)]
	    |   ('+' >> factor)
	    |   function_expression
	    |   operator_expression
	    ;
	 
	 // power operator, next precedence, operates to the right
	 pow_term =
	    factor
	    >> *(  ('^' >> pow_term)[invoke_binary<ElementType, binary_power<ElementType> >(self.eval)]
		   )
	    ;
	 
	 term =
	    pow_term
	    >> *(   ('*' >> pow_term)[invoke_binary<ElementType, 
				      binary_multiplication<ElementType> >(self.eval)]
                    |   ('/' >> pow_term)[invoke_binary<ElementType, 
					  binary_division<ElementType> >(self.eval)]
                    )
	    ;
	 
	 expression =
	    term
	    >> *(  ('+' >> term)[invoke_binary<ElementType, 
				 binary_addition<ElementType> >(self.eval)]
		   |   ('-' >> term)[invoke_binary<ElementType, 
				     binary_subtraction<ElementType> >(self.eval)]
		   )
	    >> !end_p     // skip trailing whitespace
	    ;
      }
      
      rule<ScannerT> expression, term, factor, real, imag, operator_literal, unary_function,
	 binary_function, bracket_expr, quantumnumber, prod_expression, sq_bracket_expr, 
	 operator_expression, operator_bracket_sq, operator_sq_bracket, operator_bracket, operator_sq,
	 parameter, named_parameter, parameter_list, expression_string, 
	 sum_unit_expression, sum_kink_expression, sum_k_expression,
	    identifier, pow_term, commutator_bracket, num_cells, function_expression,
	 string_expression, prod_unit_expression, prod_unit_r_expression,
	 sum_string_inner_expression, sum_string_dot_expression;

      rule<ScannerT> const& start() const { return expression; }
   };
   
   std::stack<ElementType>& eval;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   IdentifierStackType& IdentifierStack;
   FunctionStackType& FunctionStack;
   ParameterStackType& ParameterStack;
   ArgumentType& Arguments;
   InfiniteLattice const& Lattice;
};

// global variables (static members of InfiniteLatticeParser)
constants InfiniteLatticeParser::constants_p;
unary_funcs<InfiniteLatticeParser::ElementType> InfiniteLatticeParser::unary_funcs_p;
binary_funcs<InfiniteLatticeParser::ElementType> InfiniteLatticeParser::binary_funcs_p;

InfiniteMPOElement
ParseInfiniteOperator(InfiniteLattice const& Lattice, std::string const& Str,
		      Function::ArgumentList const& Args)
{
   typedef InfiniteLatticeParser::ElementType ElementType;

   InfiniteLatticeParser::ElemStackType       ElemStack;
   InfiniteLatticeParser::UnaryFuncStackType  UnaryFuncStack;
   InfiniteLatticeParser::BinaryFuncStackType BinaryFuncStack;
   InfiniteLatticeParser::IdentifierStackType IdentStack;
   InfiniteLatticeParser::ParameterStackType  ParamStack;
   InfiniteLatticeParser::FunctionStackType   FunctionStack;
   InfiniteLatticeParser::ArgumentType        Arguments;

   CheckParentheses(Str.begin(), Str.end());

   for (Function::ArgumentList::const_iterator I = Args.begin(); I != Args.end(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   // Put the lattice args into Arguments - this won't override existing values
   for (InfiniteLattice::const_argument_iterator I = Lattice.begin_arg(); I != Lattice.end_arg(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   char const* beg = Str.c_str();
   char const* end = beg + Str.size();

   InfiniteLatticeParser Parser(ElemStack, UnaryFuncStack, BinaryFuncStack, IdentStack, 
				FunctionStack, ParamStack, 
				Arguments, Lattice);
   try
   {
      parse_info<> info = parse(beg, Parser, space_p);
      if (!info.full)
	 throw ParserError::AtRange("Failed to parse an expression", info.stop, end);
   }
   catch (ParserError const& p)
   {
      throw ParserError::Finalize(p, "While parsing an infinite operator:", beg, end);
   }
   catch (std::exception const& p)
   {
      throw ParserError::Finalize(p, "While parsing an infinite operator:", beg, end);
   }
   catch (...)
   {
      throw;
   }

   CHECK(UnaryFuncStack.empty());
   CHECK(BinaryFuncStack.empty());
   CHECK(IdentStack.empty());
   CHECK(ParamStack.empty());
   CHECK(!ElemStack.empty());
   ElementType Result = ElemStack.top();
   ElemStack.pop();
   CHECK(ElemStack.empty())(ElemStack.size());

   return Result;
}

std::pair<InfiniteMPOElement, InfiniteLattice>
ParseInfiniteOperatorAndLattice(std::string const& Str)
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

   InfiniteMPOElement Op = ParseInfiniteOperator(*Lattice, Expr);
   return std::make_pair(Op, *Lattice);
}

ProductMPO
ParseProductOperator(InfiniteLattice const& Lattice, std::string const& Str,
		     Function::ArgumentList const& Args)
{
   InfiniteMPO Op = ParseInfiniteOperator(Lattice, Str, Args);
   return Op.as_product_mpo();
}

std::pair<ProductMPO, InfiniteLattice>
ParseProductOperatorAndLattice(std::string const& Str)
{
   std::pair<InfiniteMPO, InfiniteLattice>
      p = ParseInfiniteOperatorAndLattice(Str);
   return std::make_pair(p.first.as_product_mpo(), p.second);
}

TriangularMPO
ParseTriangularOperator(InfiniteLattice const& Lattice, std::string const& Str,
			Function::ArgumentList const& Args)
{
   InfiniteMPO Op = ParseInfiniteOperator(Lattice, Str, Args);
   return Op.as_triangular_mpo();
}

std::pair<TriangularMPO, InfiniteLattice>
ParseTriangularOperatorAndLattice(std::string const& Str)
{
   std::pair<InfiniteMPO, InfiniteLattice>
      p = ParseInfiniteOperatorAndLattice(Str);
   return std::make_pair(p.first.as_triangular_mpo(), p.second);
}

std::complex<double>
ParseInfiniteNumber(InfiniteLattice const& Lattice, std::string const& Str,
		    Function::ArgumentList const& Args)
{
   InfiniteMPO Op = ParseInfiniteOperator(Lattice, Str, Args);
   return Op.as_complex();
}
