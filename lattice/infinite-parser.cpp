// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/infinite-parser.cpp
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

   void operator()(char const* start, char const* end) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();
      if (!Lattice.operator_exists(OpName))
      {
         throw ParserError::AtPosition("Operator does not exist or wrong type: "
                                       + ColorQuote(OpName), start);
      }
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
                 //                 Function::ArgumentList const& Args_,
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
   //   Function::ArgumentList const& Args;
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
                  std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }
      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End), Args);
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(prod_unit_left_to_right(Op.MPO(), Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_prod_unit_r
{
   push_prod_unit_r(InfiniteLattice const& Lattice_,
                    std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End), Args);
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(prod_unit_right_to_left(Op.MPO(), Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_trans_right
{
   push_trans_right(InfiniteLattice const& Lattice_,
                    std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      std::vector<BasisList> LocalBasisList;
      for (int i = 0; i < Lattice.GetUnitCell().size(); ++i)
      {
         LocalBasisList.push_back(Lattice.GetUnitCell().LocalBasis(i).Basis());
      }
      eval.push(pow(translate_right(LocalBasisList), Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_string
{
   push_string(InfiniteLattice const& Lattice_, std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      BasicFiniteMPO Op = ParseStringOperator(*Lattice.GetUnitCell().GetSiteList(),
                                         std::string(Start, End), Lattice.GetUnitCell().size() /* , Args */);
      // TODO: Add Args to ParseStringOperator()
      eval.push(prod_unit_left_to_right(Op, Lattice.GetUnitCell().size()));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_sum_unit
{
   push_sum_unit(InfiniteLattice const& Lattice_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }

      try
      {
         UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End), Args, &Lattice);
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
   Function::ArgumentList const& Args;
};

struct push_sum_k
{
   push_sum_k(InfiniteLattice const& Lattice_,
              std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      std::complex<double> k = boost::get<std::complex<double> >(eval.top());
      eval.pop();
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }

      try
      {
         UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, End), Args);
         //Op.ExtendToCoverUnitCell(Sites);
         eval.push(sum_k(k, Op, Sites));
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
   Function::ArgumentList const& Args;
};

struct push_sum_kink
{
   push_sum_kink(InfiniteLattice const& Lattice_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }

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
         throw ParserError::AtRange(ColorQuote("sum_kink") + " expects two parameters", Start, End);
      }

      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Kink = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma), Args);
      ++Comma; // skip over the comma
      UnitCellMPO Op = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Comma, End), Args);
      Op.ExtendToCoverUnitCell(Sites);
      eval.push(sum_kink(Kink, Op, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_sum_string_inner
{
   push_sum_string_inner(InfiniteLattice const& Lattice_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }

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
         throw ParserError::AtRange(ColorQuote("sum_string_inner")
                                    + " expects three parameters, only found one", Start, End);
      }

      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op1 = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma), Args);
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
         throw ParserError::AtRange(ColorQuote("sum_string_inner")
                                    + " expects three parameters, only found two", Start, End);
      }
      UnitCellMPO String = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma), Args);
      ++Comma; // skip over the comma
      UnitCellMPO Op2 = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Comma, End), Args);
      eval.push(sum_string_inner(Op1, String, Op2, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_sum_string_dot
{
   push_sum_string_dot(InfiniteLattice const& Lattice_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }
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
         throw ParserError::AtRange(ColorQuote("sum_string_dot")
                                    + " expects three parameters, only found one", Start, End);
      }

      DEBUG_TRACE("Parsing UnitCellMPO")(std::string(Start,End));
      UnitCellMPO Op1 = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma), Args);
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
         throw ParserError::AtRange(ColorQuote("sum_string_dot")
                                    + " expects three parameters, only found two", Start, End);
      }
      UnitCellMPO String = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Start, Comma), Args);
      ++Comma; // skip over the comma
      UnitCellMPO Op2 = ParseUnitCellOperator(Lattice.GetUnitCell(), 0, std::string(Comma, End), Args);
      eval.push(sum_string_dot(Op1, String, Op2, Sites));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct push_sum_partial
{
   push_sum_partial(InfiniteLattice const& Lattice_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }

      try
      {
	        BasicTriangularMPO Op = ParseTriangularOperator(Lattice, std::string(Start, End), Args);
          eval.push(sum_partial(Op, Lattice.GetUnitCell()["I"], Sites));
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
   Function::ArgumentList const& Args;
};

struct push_sum_partial3
{
   push_sum_partial3(InfiniteLattice const& Lattice_, std::stack<std::string>& IdentifierStack_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), IdentifierStack(IdentifierStack_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      int Sites = pop_int(eval);
      if (Sites % Lattice.GetUnitCell().size() != 0)
      {
         throw ParserError::AtRange("Number  of sites must be a multiple of the lattice unit cell size",
                                    Start, End);
      }

      std::string OperatorStr = IdentifierStack.top();
      IdentifierStack.pop();

      try
      {
	        BasicTriangularMPO Op = ParseTriangularOperator(Lattice, OperatorStr, Args);
          eval.push(sum_partial(Op, ParseUnitCellOperator(Lattice.GetUnitCell(), Sites, std::string(Start, End), Args), Sites));
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
   std::stack<std::string>& IdentifierStack;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};


struct push_coarse_grain
{
   push_coarse_grain(InfiniteLattice const& Lattice_,
                 std::stack<ElementType>& eval_, Function::ArgumentList const& Args_)
      : Lattice(Lattice_), eval(eval_), Args(Args_) {}

   void operator()(char const* Start, char const* End) const
   {
      ElementType Op = eval.top();
      eval.pop();

      int Sites = pop_int(eval);
      if (Sites <= 0)
         throw ParserError::AtRange(ColorHighlight("coarse_grain:")
                                    + " canot coarse-grain to zero or negative size!",
                                    Start, End);

      BasicTriangularMPO TOp = boost::get<BasicTriangularMPO>(Op);
      TOp = coarse_grain(TOp, Sites);
      eval.push(ElementType(TOp));
   }

   InfiniteLattice const& Lattice;
   std::stack<ElementType>& eval;
   Function::ArgumentList const& Args;
};

struct init_num_param_stack
{
   init_num_param_stack(std::stack<int>& Stack_) : Stack(Stack_) {}

   void operator()(char const* Start, char const* End) const
   {
      Stack.push(1);
   }

   std::stack<int>& Stack;
};


struct increment_num_param_stack
{
   increment_num_param_stack(std::stack<int>& Stack_) : Stack(Stack_) {}

   void operator()(char const* Start, char const* End) const
   {
      CHECK(!Stack.empty());
      ++Stack.top();
   }

   std::stack<int>& Stack;
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
   typedef std::stack<int>                     NumParameterStackType;
   typedef std::stack<Function::ParameterList> ParameterStackType;
   typedef symbols<complex>                    ArgumentType;
   typedef Function::ArgumentList              RawArgumentType;

   static constants                  constants_p;
   static unary_funcs_mpo<ElementType>  unary_funcs_p;
   static binary_funcs<ElementType> binary_funcs_p;

   InfiniteLatticeParser(ElemStackType& eval_,
                  UnaryFuncStackType& func_stack_,
                  BinaryFuncStackType& bin_func_stack_,
                  IdentifierStackType& IdentifierStack_,
                  FunctionStackType& Functions_,
                  NumParameterStackType& NumParameterStack_,
                  ParameterStackType& Parameters_,
                  ArgumentType& Arguments_,
                         InfiniteLattice const& Lattice_,
                         RawArgumentType const& Args_)
      : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_),
      IdentifierStack(IdentifierStack_), FunctionStack(Functions_),
      NumParameterStack(NumParameterStack_),
      ParameterStack(Parameters_), Arguments(Arguments_),
        Lattice(Lattice_), Args(Args_)
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

         // and similarly for a filename; pretend its an identifier.
         filename = lexeme_d[*(anychar_p - chset<>(","))]
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

         bracket_expression = lexeme_d[+((anychar_p - chset<>("()"))
                                        | (ch_p('(') >> bracket_expression >> ch_p(')')))];

         expression_string = lexeme_d[+((anychar_p - chset<>("(),"))
                                        | (ch_p('(') >> bracket_expression >> ch_p(')')))];

         num_cells = (eps_p((str_p("cells") | str_p("sites")) >> '=')
                      >> ((str_p("cells") >> '=' >> expression >> ',')
                          [scale_cells_to_sites(self.Lattice, self.eval)]
                          | (str_p("sites") >> '=' >> expression >> ',')))
            | eps_p[push_number_of_sites(self.Lattice, self.eval)];

         num_cells_no_comma = (eps_p((str_p("cells") | str_p("sites")) >> '=')
                               >> ((str_p("cells") >> '=' >> expression)
                                   [scale_cells_to_sites(self.Lattice, self.eval)]
                                   | (str_p("sites") >> '=' >> expression)))
            | eps_p[push_number_of_sites(self.Lattice, self.eval)];

         // ProductMPO expressions

         string_expression = str_p("string")
            >> '('
            >> expression_string[push_string(self.Lattice, self.eval, self.Args)]
            >> ')';

         prod_unit_expression = str_p("prod_unit")
            >> '('
            >> num_cells
            >> expression_string[push_prod_unit(self.Lattice, self.eval, self.Args)]
            >> ')';

         prod_unit_r_expression = str_p("prod_unit_r")
            >> '('
            >> num_cells
            >> expression_string[push_prod_unit_r(self.Lattice, self.eval, self.Args)]
            >> ')';

         trans_right_expression = str_p("trans_right")
            >> '('
            >> num_cells_no_comma[push_trans_right(self.Lattice, self.eval, self.Args)]
            >> ')';

         // BasicTriangularMPO expressions

         sum_unit_expression = str_p("sum_unit")
            >> '('
            >> num_cells
            >> expression_string[push_sum_unit(self.Lattice, self.eval, self.Args)]
            >> ')';

         sum_kink_expression = str_p("sum_kink")
            >> '('
            >> num_cells
            >> bracket_expression[push_sum_kink(self.Lattice, self.eval, self.Args)]
            >> ')';

         sum_k_expression = str_p("sum_k")
            >> '('
            >> num_cells
            >> expression >> ','
            >> bracket_expression[push_sum_k(self.Lattice, self.eval, self.Args)]
            >> ')';

         sum_string_inner_expression = str_p("sum_string_inner")
            >> '('
            >> num_cells
            >> bracket_expression[push_sum_string_inner(self.Lattice, self.eval, self.Args)]
            >> ')';

         sum_string_dot_expression = str_p("sum_string_dot")
            >> '('
            >> num_cells
            >> bracket_expression[push_sum_string_dot(self.Lattice, self.eval, self.Args)]
            >> ')';

         sum_partial_expression = str_p("sum_partial")
            >> '('
            >> num_cells >>
              ((eps_p(expression_string >> ')') >> expression_string[push_sum_partial(self.Lattice, self.eval, self.Args)] >> ')')
              |
              (expression_string[push_identifier(self.IdentifierStack)] >> ',' >>
               expression_string[push_sum_partial3(self.Lattice, self.IdentifierStack, self.eval, self.Args)] >> ')'));

         coarse_grain_expression = str_p("coarse_grain")
            >> '('
            >> (num_cells
                >> expression >> ')')[push_coarse_grain(self.Lattice, self.eval, self.Args)];

         filegrid_expression = str_p("filegrid")
            >> '('
            >> filename >> ','
            >> expression[init_num_param_stack(self.NumParameterStack)]
            >> (*(',' >> expression[increment_num_param_stack(self.NumParameterStack)]) >> ')')
               [eval_filegrid<ElementType, InfiniteLattice>(self.Lattice, self.IdentifierStack, self.NumParameterStack, self.eval, self.Args)];

         function_expression = eps_p(identifier >> '{')
            >> identifier[push_function(self.FunctionStack, self.ParameterStack)]
            >> parameter_list[eval_function(self.Lattice,
                                            self.FunctionStack,
                                            self.ParameterStack,
                                            //                                            self.Args,
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
            |   string_expression
            |   prod_unit_expression
            |   prod_unit_r_expression
            |   trans_right_expression
            |   sum_unit_expression
            |   sum_kink_expression
            |   sum_k_expression
            |   sum_partial_expression
            |   coarse_grain_expression
            |   sum_string_inner_expression
            |   sum_string_dot_expression
            |   filegrid_expression
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
                    |   ('%' >> pow_term)[invoke_binary<ElementType,
                                          binary_modulus<ElementType> >(self.eval)]
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
         parameter, named_parameter, parameter_list, expression_string, bracket_expression,
         sum_unit_expression, sum_kink_expression, sum_k_expression, filename,
         identifier, pow_term, commutator_bracket, num_cells, num_cells_no_comma, function_expression,
         string_expression, prod_unit_expression, prod_unit_r_expression, trans_right_expression,
         sum_string_inner_expression, sum_string_dot_expression, sum_partial_expression, coarse_grain_expression,
         filegrid_expression;

      rule<ScannerT> const& start() const { return expression; }
   };

   std::stack<ElementType>& eval;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   IdentifierStackType& IdentifierStack;
   NumParameterStackType& NumParameterStack;
   FunctionStackType& FunctionStack;
   ParameterStackType& ParameterStack;
   ArgumentType& Arguments;
   InfiniteLattice const& Lattice;
   RawArgumentType const& Args;
};

// global variables (static members of InfiniteLatticeParser)
constants InfiniteLatticeParser::constants_p;
unary_funcs_mpo<InfiniteLatticeParser::ElementType> InfiniteLatticeParser::unary_funcs_p;
binary_funcs<InfiniteLatticeParser::ElementType> InfiniteLatticeParser::binary_funcs_p;

InfiniteMPOElement
ParseInfiniteOperator(InfiniteLattice const& Lattice, std::string const& Str,
                      Function::ArgumentList const& ArgsOld)
{
   typedef InfiniteLatticeParser::ElementType ElementType;

   InfiniteLatticeParser::ElemStackType       ElemStack;
   InfiniteLatticeParser::UnaryFuncStackType  UnaryFuncStack;
   InfiniteLatticeParser::BinaryFuncStackType BinaryFuncStack;
   InfiniteLatticeParser::IdentifierStackType IdentStack;
   InfiniteLatticeParser::NumParameterStackType NumParameterStack;
   InfiniteLatticeParser::ParameterStackType  ParamStack;
   InfiniteLatticeParser::FunctionStackType   FunctionStack;
   InfiniteLatticeParser::ArgumentType        Arguments;

   Function::ArgumentList Args = ArgsOld;

   CheckParentheses(Str.begin(), Str.end());

   // Add the lattice constants to Args, but don't overwrite any existing values
   for (InfiniteLattice::const_argument_iterator I = Lattice.begin_arg(); I != Lattice.end_arg(); ++I)
   {
      if (Args.find(I->first) == Args.end())
      {
         Args.insert(*I);
      }
   }

   // now make a symbols version of all of the arge

   for (Function::ArgumentList::const_iterator I = Args.begin(); I != Args.end(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   char const* beg = Str.c_str();
   char const* end = beg + Str.size();

   InfiniteLatticeParser Parser(ElemStack, UnaryFuncStack, BinaryFuncStack, IdentStack,
                                FunctionStack, NumParameterStack, ParamStack,
                                Arguments, Lattice, Args);
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

std::pair<std::string, InfiniteLattice>
ParseOperatorStringAndLattice(std::string const& Str)
{
   std::string::const_iterator Delim = std::find(Str.begin(), Str.end(), ':');
   if (Delim == Str.end())
   {
      throw ParserError("expression of the form \"lattice:expression\" expected");
   }

   std::string LatticeFile = std::string(Str.begin(), Delim);
   boost::trim(LatticeFile);
   pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

   ++Delim;
   std::string Expr(Delim, Str.end());

   return {Expr, *Lattice};
}

std::pair<InfiniteMPOElement, InfiniteLattice>
ParseInfiniteOperatorAndLattice(std::string const& Str)
{
   std::pair<std::string, InfiniteLattice> OpLattice = ParseOperatorStringAndLattice(Str);
   InfiniteMPOElement Op = ParseInfiniteOperator(OpLattice.second, OpLattice.first);
   return {Op, OpLattice.second};
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

BasicTriangularMPO
ParseTriangularOperator(InfiniteLattice const& Lattice, std::string const& Str,
                        Function::ArgumentList const& Args)
{
   InfiniteMPO Op = ParseInfiniteOperator(Lattice, Str, Args);
   return Op.as_basic_triangular_mpo();
}

std::pair<BasicTriangularMPO, InfiniteLattice>
ParseTriangularOperatorAndLattice(std::string const& Str)
{
   std::pair<InfiniteMPO, InfiniteLattice>
      p = ParseInfiniteOperatorAndLattice(Str);
   return std::make_pair(p.first.as_basic_triangular_mpo(), p.second);
}

std::complex<double>
ParseInfiniteNumber(InfiniteLattice const& Lattice, std::string const& Str,
                    Function::ArgumentList const& Args)
{
   InfiniteMPO Op = ParseInfiniteOperator(Lattice, Str, Args);
   return Op.as_complex();
}
