// -*- C++ -*- $Id: operator-parser.cpp 1326 2014-03-04 10:56:37Z ianmcc $

#include "siteoperator-parser.h"
#include "parser/parser.h"

using namespace Parser;

template <typename element_type>
struct push_operator
{
   push_operator(LatticeSite const& Site_, std::stack<element_type>& eval_)
      : Site(Site_), eval(eval_) {}

   void operator()(char const* str, char const* end) const
   {
      std::string OpName(str, end);

      CHECK(Site.operator_exists(OpName))("Operator does not exist at the lattice site")(OpName);
      eval.push(element_type(Site[OpName]));
   }

   LatticeSite const& Site;
   std::stack<element_type >& eval;
};

struct SiteOperatorParser : public grammar<SiteOperatorParser>
{
   typedef boost::variant<complex, SiteOperator> element_type;
   typedef boost::function<element_type(element_type)> unary_func_type;
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;

   typedef std::stack<element_type>     ElementStackType;
   typedef std::stack<unary_func_type>  UnaryFuncStackType;
   typedef std::stack<binary_func_type> BinaryFuncStackType;

   static constants                  constants_p;
   static unary_funcs<element_type>  unary_funcs_p;
   static binary_funcs<element_type> binary_funcs_p;

   SiteOperatorParser(ElementStackType& eval_, 
		      UnaryFuncStackType& func_stack_,
		      BinaryFuncStackType& bin_func_stack_,
		      LatticeSite const& Site_)
      : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_), Site(Site_) {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(SiteOperatorParser const& self)
      {
	 real = ureal_p[push_real<element_type>(self.eval)];
	 
	 imag = lexeme_d[ureal_p >> chset<>("iIjJ")][push_imag<element_type>(self.eval)];
	 
	 identifier = lexeme_d[alpha_p >> *(alnum_p | '_')];
	 
	 bracket_expr = '(' >> *((anychar_p - chset<>("()")) | bracket_expr) >> ')';
	 
	 //lattice_operator = (identifier >> !bracket_expr)[push_operator(self.Site, self.eval)];
	 lattice_operator = (identifier)[push_operator<element_type>(self.Site, self.eval)];
	 
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
	 binary_function,
	 bracket_expr, lattice_operator, identifier, pow_term, commutator_bracket;
      rule<ScannerT> const&
      start() const { return expression; }
   };
   
   std::stack<element_type>& eval;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   LatticeSite const& Site;
};


// global variables (static members of SiteOperatorParser)
constants SiteOperatorParser::constants_p;
unary_funcs<SiteOperatorParser::element_type> SiteOperatorParser::unary_funcs_p;
binary_funcs<SiteOperatorParser::element_type> SiteOperatorParser::binary_funcs_p;

SiteOperator
ParseSiteOperator(LatticeSite const& Site, std::string const& Str)
{
   typedef SiteOperatorParser::element_type element_type;

   SiteOperatorParser::ElementStackType ElementStack;
   SiteOperatorParser::UnaryFuncStackType UnaryFuncStack;
   SiteOperatorParser::BinaryFuncStackType BinaryFuncStack;

   SiteOperatorParser Parser(ElementStack, UnaryFuncStack, BinaryFuncStack, Site);

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

   SiteOperator* Op = boost::get<SiteOperator>(&Result);
   if (Op)
   {
      return *Op;
   }
   // else, we also handle the case where the operator is a number
   complex x = boost::get<complex>(Result);
   return x*Site["I"];
}
