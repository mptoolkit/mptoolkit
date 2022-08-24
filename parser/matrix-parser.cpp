// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// parser/matrix-parser.cpp
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

#include "parser/matrix-parser.h"
#include "parser/parser.h"
#include "tensor/tensor_exponential.h"

namespace Parser
{

template <>
std::string name_of<MatrixOperator>(MatrixOperator const&)
{
   return "matrix operator";
}

typedef boost::variant<complex, MatrixOperator> MatrixElementType;

template <>
struct ElementSq<MatrixElementType> : boost::static_visitor<MatrixElementType>
{
   MatrixElementType operator()(complex c) const
   {
      return complex(c.real()*c.real() + c.imag()*c.imag());
   }

   MatrixElementType operator()(MatrixOperator const& x) const
   {
      return inner_prod(x,x);
   }
};

} // namespace Parser

using namespace Parser;

struct MatrixParser : public grammar<MatrixParser>
{
   typedef boost::variant<complex, MatrixOperator> ElementType;
   typedef boost::function<ElementType(ElementType)> unary_func_type;
   typedef boost::function<ElementType(ElementType, ElementType)> binary_func_type;

   typedef std::stack<ElementType>            ElemStackType;
   typedef std::stack<unary_func_type>         UnaryFuncStackType;
   typedef std::stack<binary_func_type>        BinaryFuncStackType;
   typedef std::stack<std::string>             IdentifierStackType;
   typedef std::stack<std::string>             FunctionStackType;
   typedef std::stack<Function::ParameterList> ParameterStackType;
   typedef symbols<complex>                    ArgumentType;
   typedef symbols<MatrixOperator>             MatrixArgumentType;

   static constants                  constants_p;
   static unary_funcs<ElementType>  unary_funcs_p;
   static binary_funcs<ElementType> binary_funcs_p;

   MatrixParser(ElemStackType& eval_,
                UnaryFuncStackType& func_stack_,
                BinaryFuncStackType& bin_func_stack_,
                IdentifierStackType& IdentifierStack_,
                FunctionStackType& Functions_,
                ParameterStackType& Parameters_,
                ArgumentType& Arguments_,
                MatrixArgumentType& MatrixArguments_)
      : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_),
        IdentifierStack(IdentifierStack_), FunctionStack(Functions_),
        ParameterStack(Parameters_), Arguments(Arguments_), MatrixArguments(MatrixArguments_)
   {}

   template <typename ScannerT>
   struct definition
   {
      definition(MatrixParser const& self)
      {
         real_t = ureal_p;
         real = real_t[push_value<ElementType>(self.eval)];

         imag_t = lexeme_d[ureal_p >> chset<>("iIjJ")];
         imag = imag_t[push_imag<ElementType>(self.eval)];

         identifier_t = lexeme_d[alpha_p >> *(alnum_p | '_')];
         identifier = identifier_t;
         //[push_identifier(self.identifier_stack)];

         // We re-use the identifier_stack for quantum numbers, and rely on the grammar rules to
         // avoid chaos!
         quantumnumber_t = lexeme_d[*(anychar_p - chset<>("()"))];

         quantumnumber = quantumnumber_t[push_identifier(self.IdentifierStack)];

         named_parameter_t = identifier_t >> '=' >> expression_t;

         named_parameter = eps_p(identifier >> '=')
            >> identifier[push_identifier(self.IdentifierStack)]
            >> '='
            >> expression[push_named_parameter<ElementType>(self.eval,
                                                             self.IdentifierStack,
                                                             self.ParameterStack)];

         parameter_t = expression_t;
         parameter = expression[push_parameter<ElementType>(self.eval, self.ParameterStack)];

         // parameter_list is a comma-separated list of parameters, may be empty
         // at least one parameter
         parameter_list_t = '{' >> (!((named_parameter_t | parameter_t) % ',')) >> '}';
         parameter_list = '{' >> (!((named_parameter | parameter) % ',')) >> '}';

         prod_expression_t = (str_p("prod") >> '(' >> expression_t >> ','
                              >> expression_t >> ',' >> quantumnumber_t >> ')');

         prod_expression = (str_p("prod") >> '(' >> expression >> ','
                            >> expression >> ',' >> quantumnumber >> ')')
            [push_prod<ElementType>(self.IdentifierStack, self.eval)];

         bracket_expr_t = '(' >> expression_t >> ')';

         bracket_expr = '(' >> expression >> ')';

         sq_bracket_expr_t = '[' >> expression_t >> ']';

         sq_bracket_expr = '[' >> expression >> ']';

         unary_function_t = unary_funcs_p >> '(' >> expression_t >> ')';
         unary_function_t = unary_funcs_p >> '(' >> expression_t >> ')';

         unary_function =
            eps_p(unary_funcs_p >> '(')
            >>  unary_funcs_p[push_unary<ElementType>(self.func_stack)]
            >>  ('(' >> expression >> ')')[eval_unary<ElementType>(self.func_stack, self.eval)];

         binary_function_t = binary_funcs_p >> '(' >> expression_t >> ',' >> expression_t >> ')';

         binary_function =
            eps_p(self.binary_funcs_p >> '(')
            >>  self.binary_funcs_p[push_binary<ElementType>(self.bin_func_stack)]
            >>  ('(' >> expression >> ','  >> expression >> ')')
            [eval_binary<ElementType>(self.bin_func_stack, self.eval)];

         commutator_bracket_t = ch_p('[') >> expression_t >> ']';

         commutator_bracket =
            ('[' >> expression >> ',' >> expression >> ']')[invoke_binary<ElementType,
                                                            binary_commutator<ElementType> >(self.eval)];


         factor_t =
            imag_t
            |   real_t
            |   unary_function_t
            |   binary_function_t
            |   keyword_d[constants_p[eps_p]]
            |   keyword_d[self.Arguments]
            |   keyword_d[self.MatrixArguments]
            |   prod_expression_t
            |   commutator_bracket_t
            |   '(' >> expression_t >> ')'
            |   ('-' >> factor_t)
            |   ('+' >> factor_t)
            ;

         factor =
            imag
            |   real
            |   unary_function
            |   binary_function
            |   keyword_d[constants_p[push_value<ElementType>(self.eval)]]
            |   keyword_d[self.Arguments[push_value<ElementType>(self.eval)]]
            |   keyword_d[self.MatrixArguments[push_value<ElementType>(self.eval)]]
            |   prod_expression
            |   commutator_bracket
            |   '(' >> expression >> ')'
            |   ('-' >> factor)[do_negate<ElementType>(self.eval)]
            |   ('+' >> factor)
            ;

         // power operator, next precedence, operates to the right
         pow_term_t =
            factor_t
            >> *(('^' >> pow_term_t));

         pow_term =
            factor
            >> *(  ('^' >> pow_term)[invoke_binary<ElementType, binary_power<ElementType> >(self.eval)]
                   )
            ;

         term_t =
            pow_term_t
            >> *(   ('*' >> pow_term_t)
                    |   ('/' >> pow_term_t)
                    |   ('%' >> pow_term_t)
                    )
            ;
         term =
            pow_term
            >> *(   ('*' >> pow_term)[invoke_binary<ElementType,
                                      binary_multiplication<ElementType> >(self.eval)]
                    |   ('/' >> pow_term)[invoke_binary<ElementType,
                                          binary_division<ElementType> >(self.eval)]
                    |   ('%' >> pow_term)[invoke_binary<ElementType,
                                          binary_modulus<ElementType> >(self.eval)]
                    )
            ;

         expression_t =
            term_t
            >> *(  ('+' >> term_t)
                   |   ('-' >> term_t)
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

      rule<ScannerT> real, imag, identifier, quantumnumber,
         parameter, named_parameter, parameter_list,
         prod_expression,
         bracket_expr, sq_bracket_expr,
         string_expression,
         unary_function, binary_function,
         commutator_bracket,
         factor, pow_term, term, expression;

      rule<ScannerT> real_t, imag_t, identifier_t, quantumnumber_t,
         parameter_t, named_parameter_t, parameter_list_t,
         prod_expression_t,
         bracket_expr_t, sq_bracket_expr_t,
         string_expression_t,
         unary_function_t, binary_function_t,
         commutator_bracket_t,
         factor_t, pow_term_t, term_t, expression_t;

      rule<ScannerT> const& start() const { return expression; }
   };

   std::stack<ElementType>& eval;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   IdentifierStackType& IdentifierStack;
   FunctionStackType& FunctionStack;
   ParameterStackType& ParameterStack;
   ArgumentType& Arguments;
   MatrixArgumentType& MatrixArguments;
};


// global variables (static members of MatrixParser)
constants MatrixParser::constants_p;
unary_funcs<MatrixParser::ElementType> MatrixParser::unary_funcs_p;
binary_funcs<MatrixParser::ElementType> MatrixParser::binary_funcs_p;

MatrixOperator
ParseMatrixOperator(std::string const& Str,
                    Function::ArgumentList const& Args,
                    std::map<std::string, MatrixOperator> const& Matrices)
{
   typedef MatrixParser::ElementType ElementType;

   MatrixParser::ElemStackType       ElemStack;
   MatrixParser::UnaryFuncStackType  UnaryFuncStack;
   MatrixParser::BinaryFuncStackType BinaryFuncStack;
   MatrixParser::IdentifierStackType IdentStack;
   MatrixParser::ParameterStackType  ParamStack;
   MatrixParser::FunctionStackType   FunctionStack;
   MatrixParser::ArgumentType        Arguments;
   MatrixParser::MatrixArgumentType  MatrixArguments;

   CheckParentheses(Str.begin(), Str.end());

   for (Function::ArgumentList::const_iterator I = Args.begin(); I != Args.end(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   for (auto const& x : Matrices)
   {
      MatrixArguments.add(x.first.c_str(), x.second);
   }

   char const* beg = Str.c_str();
   char const* end = beg + Str.size();

   MatrixParser Parser(ElemStack, UnaryFuncStack, BinaryFuncStack, IdentStack,
                       FunctionStack, ParamStack,
                       Arguments, MatrixArguments);

   try
   {
      parse_info<> info = parse(beg, Parser, space_p);
      if (!info.full)
         throw ParserError::AtRange("Failed to parse an expression", info.stop, end);
   }
   catch (ParserError const& p)
   {
      throw ParserError::Finalize(p, "While parsing a matrix expression:", beg, end);
   }
   catch (std::exception const& p)
   {
      throw ParserError::Finalize(p, "While parsing a matrix expression:", beg, end);
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

   if (boost::get<MatrixOperator>(&Result))
      return boost::get<MatrixOperator>(Result);
   else
   {
      std::map<std::string, MatrixOperator>::const_iterator IPtr = Matrices.find("I");
      CHECK(IPtr != Matrices.end());
      return boost::get<std::complex<double>>(Result) * IPtr->second;
   }
}
