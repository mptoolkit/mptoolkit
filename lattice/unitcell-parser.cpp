// -*- C++ -*- $Id$

#include "unitcell-parser.h"
#include "parser/parser.h"

using namespace Parser;

typedef boost::variant<complex, FiniteMPO> element_type;

struct push_operator
{
   push_operator(UnitCell const& Cell_, 
		 std::stack<std::string>& IdentifierStack_, 
		 std::stack<element_type>& eval_)
      : Cell(Cell_), IdentifierStack(IdentifierStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      CHECK(Cell.operator_exists(OpName))("Operator does not exist in the unit cell")(OpName);
      eval.push(element_type(Cell.Operator(OpName)));
   }

   UnitCell const& Cell;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
};

struct push_local_operator
{
   push_local_operator(UnitCell const& Cell_, 
		       std::stack<std::string>& IdentifierStack_, 
		       std::stack<element_type>& eval_)
      : Cell(Cell_), IdentifierStack(IdentifierStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      complex SiteIndex = boost::get<complex>(eval.top());
      eval.pop();
      int n = int(SiteIndex.real() + 0.5);
      CHECK(norm_frob(SiteIndex - double(n)) < 1E-7)("Site index must be an integer")(SiteIndex);

      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);
      eval.push(element_type(Cell.Operator(OpName, n)));
   }

   UnitCell const& Cell;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
};


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

struct UnitCellParser : public grammar<UnitCellParser>
{
   typedef boost::variant<complex, FiniteMPO> element_type;
   typedef boost::function<element_type(element_type)> unary_func_type;
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;

   typedef std::stack<element_type>     ElementStackType;
   typedef std::stack<std::string>      IdentifierStackType;
   typedef std::stack<unary_func_type>  UnaryFuncStackType;
   typedef std::stack<binary_func_type> BinaryFuncStackType;

   static constants                  constants_p;
   static unary_funcs<element_type>  unary_funcs_p;
   static binary_funcs<element_type> binary_funcs_p;

   UnitCellParser(ElementStackType& eval_, 
		  IdentifierStackType& identifier_stack_,
		  UnaryFuncStackType& func_stack_,
		  BinaryFuncStackType& bin_func_stack_,
		  UnitCell const& Cell_)
      : eval(eval_), identifier_stack(identifier_stack_), 
	func_stack(func_stack_), bin_func_stack(bin_func_stack_), Cell(Cell_) {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(UnitCellParser const& self)
      {
	 real = ureal_p[push_real<element_type>(self.eval)];
	 
	 imag = lexeme_d[ureal_p >> chset<>("iIjJ")][push_imag<element_type>(self.eval)];
	 
	 identifier = lexeme_d[alpha_p >> *(alnum_p | '_')][push_identifier(self.identifier_stack)];
	 
	 bracket_expr = '(' >> *((anychar_p - chset<>("()")) | bracket_expr) >> ')';

	 sq_bracket_expr = '[' >> expression >> ']';
	 
	 lattice_operator = identifier
	    [push_operator(self.Cell, self.identifier_stack, self.eval)];
	 
	 local_operator = (identifier >> sq_bracket_expr)
	    [push_local_operator(self.Cell, self.identifier_stack, self.eval)];
	 
	 unary_function = 
	    eps_p(unary_funcs_p >> '(') 
	    >>  unary_funcs_p[push_unary<element_type>(self.func_stack)]
	    >>  ('(' >> expression >> ')')[eval_unary<element_type>(self.func_stack, self.eval)];
	 
	 binary_function = 
	    eps_p(binary_funcs_p >> '(') 
	    >>  binary_funcs_p[push_binary<element_type>(self.bin_func_stack)]
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
	    |   commutator_bracket
	    |   '(' >> expression >> ')'
	    |   ('-' >> factor)[do_negate<element_type>(self.eval)]
	    |   ('+' >> factor)
	    |   local_operator
	    |   lattice_operator
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
	    binary_function, bracket_expr, sq_bracket_expr, lattice_operator, 
	    local_operator, identifier, pow_term, commutator_bracket;
      rule<ScannerT> const&
      start() const { return expression; }
   };
   
   std::stack<element_type>& eval;
   IdentifierStackType& identifier_stack;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   UnitCell const& Cell;
};


// global variables (static members of UnitCellParser)
constants UnitCellParser::constants_p;
unary_funcs<UnitCellParser::element_type> UnitCellParser::unary_funcs_p;
binary_funcs<UnitCellParser::element_type> UnitCellParser::binary_funcs_p;

FiniteMPO
ParseUnitCellOperator(UnitCell const& Cell, std::string const& Str)
{
   typedef UnitCellParser::element_type element_type;

   UnitCellParser::ElementStackType ElementStack;
   UnitCellParser::IdentifierStackType IdentifierStack;
   UnitCellParser::UnaryFuncStackType UnaryFuncStack;
   UnitCellParser::BinaryFuncStackType BinaryFuncStack;

   UnitCellParser Parser(ElementStack, IdentifierStack, UnaryFuncStack, BinaryFuncStack, Cell);

   parse_info<> info = parse(Str.c_str(), Parser, space_p);
   if (!info.full)
   {
      PANIC("Operator parser failed, stopped at")(info.stop);
   }

   CHECK(UnaryFuncStack.empty());
   CHECK(BinaryFuncStack.empty());
   CHECK(!ElementStack.empty());
   element_type Result = ElementStack.top();
   ElementStack.pop();
   CHECK(ElementStack.empty());

   FiniteMPO* Op = boost::get<FiniteMPO>(&Result);
   if (Op)
   {
      return *Op;
   }
   // else, we also handle the case where the operator is a number
   complex x = boost::get<complex>(Result);
   return x*Cell.Operator("I");
}
