// -*- C++ -*- $Id: operator-parser.cpp 1326 2014-03-04 10:56:37Z ianmcc $

#include "siteoperator-parser.h"
#include "parser/parser.h"
#include "function.h"

using namespace Parser;

typedef boost::variant<SiteOperator, std::complex<double> > element_type;

namespace Parser
{

template <>
std::string name_of<SiteOperator>(SiteOperator const&)
{
   return "local operator";
}

}

struct push_site_operator
{
   push_site_operator(std::stack<element_type>& eval_, LatticeSite const& Site_)
      : eval(eval_), Site(Site_) {}

   void operator()(char const* str, char const* end) const
   {
      std::string OpName(str, end);

      CHECK(Site.operator_exists(OpName))("Operator does not exist at the lattice site")(OpName);
      eval.push(element_type(Site[OpName]));
   }

   std::stack<element_type>& eval;
   LatticeSite const& Site;
};

struct eval_function
{
   eval_function(std::stack<std::string>& FuncStack_, 
		 std::stack<Function::ParameterList>& ParamStack_,
		 std::stack<element_type>& eval_,
		 LatticeSite const& Site_)
      : FuncStack(FuncStack_), ParamStack(ParamStack_), eval(eval_), Site(Site_) {}

   void operator()(char const*, char const*) const
   {
      eval.push(Site.eval_function(FuncStack.top(), ParamStack.top()));
      FuncStack.pop();
      ParamStack.pop();
   }

   std::stack<std::string>& FuncStack;
   std::stack<Function::ParameterList>& ParamStack;
   std::stack<element_type>& eval;
   LatticeSite const& Site;
};

struct SiteOperatorParser : public grammar<SiteOperatorParser>
{
   typedef boost::function<element_type(element_type)> unary_func_type;
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;

   typedef std::stack<element_type>            ElementStackType;
   typedef std::stack<unary_func_type>         UnaryFuncStackType;
   typedef std::stack<binary_func_type>        BinaryFuncStackType;
   typedef std::stack<std::string>             IdentifierStackType;
   typedef std::stack<std::string>             FunctionStackType;
   typedef std::stack<Function::ParameterList> ParameterStackType;
   typedef symbols<complex>                    ArgumentType;

   static constants                  constants_p;
   static unary_funcs<element_type>  unary_funcs_p;
   static binary_funcs<element_type> binary_funcs_p;

   SiteOperatorParser(ElementStackType& eval_, 
		      UnaryFuncStackType& func_stack_,
		      BinaryFuncStackType& bin_func_stack_,
		      IdentifierStackType& IdentifierStack_,
		      FunctionStackType& Functions_,
		      ParameterStackType& Parameters_,
		      ArgumentType& Arguments_,
		      LatticeSite const& Site_)
      : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_), 
	IdentifierStack(IdentifierStack_), FunctionStack(Functions_),
	ParameterStack(Parameters_), 
	Arguments(Arguments_), Site(Site_) {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(SiteOperatorParser const& self)
      {
	 real = ureal_p[push_value<element_type>(self.eval)];
	 
	 imag = lexeme_d[ureal_p >> chset<>("iIjJ")][push_imag<element_type>(self.eval)];
	 
	 identifier = lexeme_d[alpha_p >> *(alnum_p | '_')];
	 
	 bracket_expr = '(' >> *((anychar_p - chset<>("()")) | bracket_expr) >> ')';
	 
	 lattice_operator = identifier[push_site_operator(self.eval, self.Site)];

	 named_parameter = eps_p(identifier >> '=')
	    >> identifier[push_identifier(self.IdentifierStack)]
	    >> '='
	    >> expression[push_named_parameter<element_type>(self.eval, 
							     self.IdentifierStack, 
							     self.ParameterStack)];

	 parameter = expression[push_parameter<element_type>(self.eval, self.ParameterStack)];

	 // parameter_list is a comma-separated list of parameters, may be empty
	 // at least one parameter
	 parameter_list = !((named_parameter | parameter) % ',');
	 
	 siteoperator_function = eps_p(identifier >> '{')
	    >> identifier[push_function(self.FunctionStack, self.ParameterStack)]
	    >> ('{' >> parameter_list >> '}')[eval_function(self.FunctionStack,
							    self.ParameterStack,
							    self.eval,
							    self.Site)];

	 unary_function = 
	    eps_p(unary_funcs_p >> '(') 
	    >>  unary_funcs_p[push_unary<element_type>(self.func_stack)]
	    >>  ('(' >> expression >> ')')[eval_unary<element_type>(self.func_stack, 
								    self.eval)];
	 
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
	    |   keyword_d[constants_p[push_value<element_type>(self.eval)]]
	    |   keyword_d[self.Arguments[push_value<element_type>(self.eval)]]
	    |   commutator_bracket
	    |   '(' >> expression >> ')'
	    |   ('-' >> factor)[do_negate<element_type>(self.eval)]
	    |   ('+' >> factor)
	    |   siteoperator_function
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
	    >> !end_p     // skip trailing whitespace
	    ;
      }
      
      rule<ScannerT> expression, term, factor, real, imag, operator_literal, unary_function,
	 binary_function, parameter, named_parameter, parameter_list, siteoperator_function,
	 bracket_expr, lattice_operator, identifier, pow_term, commutator_bracket;
      rule<ScannerT> const&
      start() const { return expression; }
   };
   
   ElementStackType& eval;
   UnaryFuncStackType& func_stack;
   BinaryFuncStackType& bin_func_stack;
   IdentifierStackType& IdentifierStack;
   FunctionStackType& FunctionStack;
   ParameterStackType& ParameterStack;
   ArgumentType& Arguments;
   LatticeSite const& Site;
};

// global variables (static members of SiteOperatorParser)
constants SiteOperatorParser::constants_p;
unary_funcs<element_type> SiteOperatorParser::unary_funcs_p;
binary_funcs<element_type> SiteOperatorParser::binary_funcs_p;

SiteElementType
ParseSiteElement(LatticeSite const& Site, 
		 std::string const& Str, 
		 Function::ArgumentList const& Args)
{
   SiteOperatorParser::ElementStackType    ElementStack;
   SiteOperatorParser::UnaryFuncStackType  UnaryFuncStack;
   SiteOperatorParser::BinaryFuncStackType BinaryFuncStack;
   SiteOperatorParser::IdentifierStackType IdentifierStack;
   SiteOperatorParser::FunctionStackType   FunctionStack;
   SiteOperatorParser::ParameterStackType  ParameterStack;
   SiteOperatorParser::ArgumentType        Arguments;

   CheckParentheses(Str.begin(), Str.end());

   // put Args into Arguments
   for (Function::ArgumentList::const_iterator I = Args.begin(); I != Args.end(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   // Put the lattice args into Arguments - this won't override existing values
   for (LatticeSite::const_argument_iterator I = Site.begin_arg(); I != Site.end_arg(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   char const* beg = Str.c_str();
   char const* end = beg + Str.size();

   SiteOperatorParser Parser(ElementStack, UnaryFuncStack, BinaryFuncStack, 
			     IdentifierStack, FunctionStack, ParameterStack, Arguments, Site);

   try
   {
      parse_info<> info = parse(beg, Parser, space_p);
      if (!info.full)
	 throw ParserError::AtRange("Failed to parse an expression", info.stop, end);
   }
   catch (ParserError const& p)
   {
      throw ParserError::Finalize(p, "While parsing a local operator:", beg, end);
   }
   catch (std::exception const& p)
   {
      throw ParserError::Finalize(p, "While parsing a local operator:", beg, end);
   }
   catch (...)
   {
      throw;
   }

   CHECK(!ElementStack.empty());
   CHECK(UnaryFuncStack.empty());
   CHECK(BinaryFuncStack.empty());
   CHECK(FunctionStack.empty());
   CHECK(ParameterStack.empty());
   element_type Result = ElementStack.top();
   ElementStack.pop();
   CHECK(ElementStack.empty());

   return Result;
}

SiteOperator
ParseSiteOperator(LatticeSite const& Site, 
		  std::string const& Str, 
		  Function::ArgumentList const& Args)
{
   SiteElementType Result = ParseSiteElement(Site, Str, Args);

   SiteOperator* Op = boost::get<SiteOperator>(&Result);
   if (Op)
   {
      return *Op;
   }
   // else, we also handle the case where the operator is a number
   complex x = boost::get<complex>(Result);
   return x*Site.identity();
}

std::complex<double>
ParseSiteNumber(LatticeSite const& Site, 
		std::string const& Str,
		Function::ArgumentList const& Args)
{
   SiteElementType Result = ParseSiteElement(Site, Str, Args);
   complex* x = boost::get<complex>(&Result);
   CHECK(x)("expression is an operator, but a c-number was expected!")(Str);
   return *x;
}
