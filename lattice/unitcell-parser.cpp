// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// lattice/unitcell-parser.cpp
//
// Copyright (C) 2014-2022 Ian McCulloch <ian@qusim.net>
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

#include "unitcell-parser.h"
#include "siteoperator-parser.h"
#include "infinite-parser.h"
#include "parser/parser.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>

InfiniteLattice const* ILattice = NULL;

namespace Parser
{
template <>
std::string name_of<UnitCellMPO>(UnitCellMPO const&)
{
   return "unit cell operator";
}

}

namespace UP
{

using namespace Parser;

typedef boost::variant<complex, UnitCellMPO> ElementType;

// Operator expressions:
// Op(cell)[site] - a local operator at a given site of a given unit cell
// Op(cell)       - a cell operator at a given cell index

// Function expressions:
// f(cell)[site]{params} - a local function
// f(cell){params}       - a UnitCell function

// Argument expressions:
// Arg[site]             - a local argument
// Arg                   - a UnitCell argument

// Op(cell)[site]
// The cell index was pushed first, in the notation Op(n)[j]
struct push_local_operator
{
   push_local_operator(UnitCell const& Cell_, int NumCells_,
                       std::stack<std::string>& IdentStack_,
                       std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const* start, char const* end) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int n = pop_int(eval);
      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Site index out of bounds")(j)(NumCells);
      if (!Cell.local_operator_exists(OpName, n))
      {
         ParserError w = ParserError::AtRange("Local operator does not exist: '"
                                              + ColorHighlight(OpName) + "'", start, end);
         if (Cell.operator_exists(OpName))
         {
            w.AddHint("Did you intend the unit cell operator "
                      + ColorHighlight(OpName + "(" +
                                      boost::lexical_cast<std::string>(n) + ")") + " ?");
         }
         throw w;
      }

      // Fetch the operator and JW string
      eval.push(ElementType(Cell.local_operator(OpName, j, n)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<ElementType >& eval;
};

// Op(cell)
// The cell number j must satisfy j >= 0 && (NumCells == 0 || j < NumCells)
struct push_cell_operator
{
   push_cell_operator(UnitCell const& Cell_, int NumCells_,
                      std::stack<std::string>& IdentStack_,
                      std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const* a, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      try
      {
         int j = pop_int(eval);
         if (!Cell.operator_exists(OpName))
            throw ParserError("Operator not found in the unit cell: '" + ColorHighlight(OpName)
                              + "'");
         eval.push(ElementType(Cell(OpName, j)));
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

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<ElementType >& eval;
};

//
// functions
//

// function f(cell)[site]{params}
struct eval_local_function
{
   eval_local_function(UnitCell const& Cell_,
                       std::stack<std::string>& FuncStack_,
                       std::stack<Function::ParameterList>& ParamStack_,
                       std::stack<ElementType>& eval_)
      : Cell(Cell_), FuncStack(FuncStack_),
        ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      TRACE(eval.size());
      int Site = pop_int(eval);
      int n = pop_int(eval);

      TRACE(n)(Site)(&Cell);

      SiteElementType Element = Cell[Site].eval_function(FuncStack.top(), ParamStack.top());
      FuncStack.pop();
      ParamStack.pop();

      // If we have a c-number, we can return it immediately
      std::complex<double>* c = boost::get<std::complex<double> >(&Element);
      if (c)
      {
         eval.push(*c);
      }
      else
      {
         SiteOperator Op = boost::get<SiteOperator>(Element);
         eval.push(Cell.map_local_operator(Op, n, Site));
      }
      TRACE(eval.size());
   }

   UnitCell const& Cell;
   std::stack<std::string>& FuncStack;
   std::stack<Function::ParameterList>& ParamStack;
   std::stack<ElementType >& eval;
};

// a c-number function f[site]{params}
struct eval_local_c_function
{
   eval_local_c_function(UnitCell const& Cell_,
                         std::stack<std::string>& FuncStack_,
                         std::stack<Function::ParameterList>& ParamStack_,
                         std::stack<ElementType>& eval_)
      : Cell(Cell_), FuncStack(FuncStack_),
        ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      int Site = pop_int(eval);
      SiteElementType Result = Cell[Site].eval_function(FuncStack.top(), ParamStack.top());

      std::complex<double>* c = boost::get<std::complex<double> >(&Result);
      CHECK(c)("Local function must return a c-number")(FuncStack.top())(Site);

      FuncStack.pop();
      ParamStack.pop();

      eval.push(*c);
   }

   UnitCell const& Cell;
   std::stack<std::string>& FuncStack;
   std::stack<Function::ParameterList>& ParamStack;
   std::stack<ElementType >& eval;
};

// f(cell){params}
struct eval_cell_function
{
   eval_cell_function(UnitCell const& Cell_,
                      int NumCells_,
                      std::stack<std::string>& FuncStack_,
                      std::stack<Function::ParameterList>& ParamStack_,
                      std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), FuncStack(FuncStack_),
        ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Cell index out of bounds")(j)(NumCells);

      // Construct the operator
      UnitCellElementType Result = Cell.eval_function(FuncStack.top(), j, ParamStack.top());

      FuncStack.pop();
      ParamStack.pop();

      eval.push(Result);
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& FuncStack;
   std::stack<Function::ParameterList>& ParamStack;
   std::stack<ElementType >& eval;
};

// a c-number function of the form f{Params}
struct eval_cell_c_function
{
   eval_cell_c_function(UnitCell const& Cell_,
                      std::stack<std::string>& FuncStack_,
                      std::stack<Function::ParameterList>& ParamStack_,
                      std::stack<ElementType>& eval_)
      : Cell(Cell_), FuncStack(FuncStack_),
        ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const* beg, char const* end) const
   {
      // Construct the operator
      UnitCellElementType Result = Cell.eval_function(FuncStack.top(), 0, ParamStack.top());

      // the result must be a c-number

      std::complex<double>* c = boost::get<std::complex<double> >(&Result);
      if (!c)
      {
         throw ParserError::AtRange("Function without cell index must return a c-number: "
                                    + ColorHighlight(FuncStack.top()),
                                    beg, end);
      }

      FuncStack.pop();
      ParamStack.pop();

      eval.push(*c);
   }

   UnitCell const& Cell;
   std::stack<std::string>& FuncStack;
   std::stack<Function::ParameterList>& ParamStack;
   std::stack<ElementType >& eval;
};

// local argument

struct push_local_argument
{
   push_local_argument(UnitCell const& Cell_,
                       std::stack<std::string>& IdentStack_,
                       std::stack<ElementType>& eval_)
      : Cell(Cell_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string ArgName = IdentStack.top();
      IdentStack.pop();

      int n = pop_int(eval);
      CHECK(n >= 0 && n < Cell.size())("Site index is out of bounds")(n)(Cell.size());

      eval.push(ElementType(Cell[n].arg(ArgName)));
   }

   UnitCell const& Cell;
   std::stack<std::string>& IdentStack;
   std::stack<ElementType >& eval;
};

//
// swapb - swap without fermion signs
//

// notation for swapb(cell1,cell2)[site1,site2]
struct push_swapb_cell_site
{
   push_swapb_cell_site(UnitCell const& Cell_, int NumCells_,
                       std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);
      int Cell2 = pop_int(eval);
      int Cell1 = pop_int(eval);

      CHECK(NumCells == 0 || (Cell1 >= 0 && Cell1 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);
      CHECK(NumCells == 0 || (Cell2 >= 0 && Cell2 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);

      eval.push(ElementType(Cell.swap_gate_no_sign(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

// notation for swapb[site1,site2]
struct push_swapb_site
{
   push_swapb_site(UnitCell const& Cell_, int NumCells_,
                  std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      CHECK(NumCells == 0 || NumCells == 1)("Cell index required");
      int Cell2 = 0;
      int Cell1 = 0;
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);

      eval.push(ElementType(Cell.swap_gate_no_sign(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

// notation for swapb(cell1,cell2)
struct push_swapb_cell
{
   push_swapb_cell(UnitCell const& Cell_, int NumCells_,
                  std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      CHECK_EQUAL(Cell.size(), 1)("Unit cell is more than one site, so a site index required here");
      int Cell2 = pop_int(eval);
      int Cell1 = pop_int(eval);
      int Site2 = 0;
      int Site1 = 0;

      CHECK(NumCells == 0 || (Cell1 >= 0 && Cell1 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);
      CHECK(NumCells == 0 || (Cell2 >= 0 && Cell2 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);

      eval.push(ElementType(Cell.swap_gate_no_sign(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

//
// swap (with fermion signs)
//

// notation for swap(cell1,cell2)[site1,site2]
struct push_swap_cell_site
{
   push_swap_cell_site(UnitCell const& Cell_, int NumCells_,
                       std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);
      int Cell2 = pop_int(eval);
      int Cell1 = pop_int(eval);

      CHECK(NumCells == 0 || (Cell1 >= 0 && Cell1 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);
      CHECK(NumCells == 0 || (Cell2 >= 0 && Cell2 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);

      PANIC("Swap() is not yet implemented");
      eval.push(ElementType(Cell.swap_gate(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

// notation for swap[site1,site2]
struct push_swap_site
{
   push_swap_site(UnitCell const& Cell_, int NumCells_,
                  std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      CHECK(NumCells == 0 || NumCells == 1)("Cell index required");
      int Cell2 = 0;
      int Cell1 = 0;
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);

      PANIC("Swap() is not yet implemented");
      eval.push(ElementType(Cell.swap_gate(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

// notation for swapb(cell1,cell2)
struct push_swap_cell
{
   push_swap_cell(UnitCell const& Cell_, int NumCells_,
                  std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      CHECK_EQUAL(Cell.size(), 1)("Unit cell is more than one site, so a site index required here");
      int Cell2 = pop_int(eval);
      int Cell1 = pop_int(eval);
      int Site2 = 0;
      int Site1 = 0;

      CHECK(NumCells == 0 || (Cell1 >= 0 && Cell1 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);
      CHECK(NumCells == 0 || (Cell2 >= 0 && Cell2 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);

      PANIC("Swap() is not yet implemented");
      eval.push(ElementType(Cell.swap_gate(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

struct push_string
{
   push_string(UnitCell const& Cell_, int NumCells_,
               std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_)
   {
      // If the number of unit cells is unspecified, use 1
      if (NumCells == 0)
         NumCells = 1;
   }

   void operator()(char const* Start, char const* End) const
   {
      BasicFiniteMPO Op = ParseStringOperator(*Cell.GetSiteList(),
                                         std::string(Start, End), NumCells*Cell.size());
      eval.push(UnitCellMPO(Cell.GetSiteList(), Op, LatticeCommute::Bosonic, 0));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

BasicFiniteMPO
fsup(int Offset1, int Offset2, BasicTriangularMPO const& Op)
{
   int const Size = Offset2 - Offset1;
   Offset1 = Offset1 % Op.size();
   GenericMPO Result(Size);
   for (int i = 0; i < Size; ++i)
   {
      Result[i] = Op[(i+Offset1)%Op.size()];
   }
   // patch up the first and last sites
   Result.front() = project_rows(Result.front(), std::set<int>({0}));
   Result.back() = project_columns(Result.back(), std::set<int>({int(Result.back().Basis2().size())-1}));
   return BasicFiniteMPO(Result);
}

BasicFiniteMPO
fsup(int Offset1, int Offset2, ProductMPO const& Op)
{
   PANIC("not yet implemented");
}

BasicFiniteMPO
fsup(int Offset1, int Offset2, InfiniteMPOElement const& Op)
{
   if (boost::get<BasicTriangularMPO>(&Op))
   {
      return fsup(Offset1, Offset2, boost::get<BasicTriangularMPO>(Op));
   }
   else if (boost::get<ProductMPO>(&Op))
   {
      return fsup(Offset1, Offset2, boost::get<ProductMPO>(Op));
   }
   else
   {
      PANIC("Unknown operator type");
   }
}

struct push_fsup
{
   push_fsup(UnitCell const& Cell_, int NumCells_,
             std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_)
   {
      // If the number of unit cells is unspecified, use 1
      if (NumCells == 0)
         NumCells = 1;
   }

   void operator()(char const* Start, char const* End) const
   {
      if (!ILattice)
      {
         PANIC("Lattice not defined.");
      }
      int Cell2 = pop_int(eval);
      int Cell1 = pop_int(eval);
      std::string Str(Start, End);
      //      TRACE(Str);
      InfiniteMPOElement Op = ParseInfiniteOperator(*ILattice, Str);

      BasicFiniteMPO NewOp = fsup(Cell1*Cell.size(), Cell2*Cell.size(), Op);
      eval.push(UnitCellMPO(Cell.GetSiteList(), NewOp, LatticeCommute::Bosonic, Cell1*Cell.size()));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType>& eval;
};

// convert the top of the expression stack from a number of cells
// to a number of sites, by multiplying by the unit cell size
struct scale_cells_to_sites
{
   scale_cells_to_sites(UnitCell const& Cell_, int NumCells_,
                        std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const* a, char const*) const
   {
      try
      {
         int Cells = pop_int(eval);
         eval.push(std::complex<double>(double(Cells*Cell.size())));
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

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};

// push the number of sites in the unit cell onto the top of the expression stack
struct push_number_of_sites
{
   push_number_of_sites(UnitCell const& Cell_, int NumCells_,
                        std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      eval.push(std::complex<double>(double(Cell.size())));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType >& eval;
};


struct push_coarse_grain
{
   push_coarse_grain(UnitCell const& Cell_, int NumCells_,
                 std::stack<ElementType>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const* Start, char const* End) const
   {
      ElementType Op = eval.top();
      eval.pop();

      int Sites = pop_int(eval);
      if (Sites <= 0)
         throw ParserError::AtRange(ColorHighlight("coarse_grain:")
                                    + " canot coarse-grain to zero or negative size!",
                                    Start, End);

      UnitCellMPO TOp = boost::get<UnitCellMPO>(Op);
      TOp = coarse_grain(TOp, Sites);
      eval.push(ElementType(TOp));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<ElementType>& eval;
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

struct push_default_randexpr
{
   push_default_randexpr(std::stack<ElementType>& eval_)
      : eval(eval_) {}

   void operator()(char const* Start, char const* End) const
   {
      eval.push(std::complex<double>(0.0));
      eval.push(std::complex<double>(1.0));
   }
   std::stack<ElementType>& eval;
};

struct push_rand_expr
{
   push_rand_expr(std::stack<int>& NumParam_, std::stack<ElementType>& eval_)
      : NumParam(NumParam_), eval(eval_) {}

   void operator()(char const* Start, char const* End) const
   {
      int n = NumParam.top();
      NumParam.pop();

      std::vector<uint32_t> Data(n);
      for (int i = 0; i < n; ++i)
      {
         Data[i] = pop_int(eval);
      }
      std::array<uint32_t, 8> Digest = SHA256::hash(Data);
      uint64_t k = (uint64_t(Digest[2]) << 32) | uint64_t(Digest[6]);
      double x = k * 5.42101086242752217e-20;

      double Max = pop_real(eval);
      double Min = pop_real(eval);
      if (Max <= Min)
      {
         throw ParserError::AtRange(ColorHighlight("rand:")
                                    + " Range is not valid, require min < max!",
                                    Start, End);
      }
      eval.push(std::complex<double>(x * (Max - Min) + Min));
   }

   std::stack<int>& NumParam;
   std::stack<ElementType>& eval;
};


} // namespace UP

using namespace UP;

struct UnitCellParser : public grammar<UnitCellParser>
{
   typedef boost::variant<complex, UnitCellMPO> ElementType;
   typedef boost::function<ElementType(ElementType)> unary_func_type;
   typedef boost::function<ElementType(ElementType, ElementType)> binary_func_type;

   typedef std::stack<ElementType>            ElemStackType;
   typedef std::stack<unary_func_type>         UnaryFuncStackType;
   typedef std::stack<binary_func_type>        BinaryFuncStackType;
   typedef std::stack<std::string>             IdentifierStackType;
   typedef std::stack<int>                     NumParameterStackType;
   typedef std::stack<std::string>             FunctionStackType;
   typedef std::stack<Function::ParameterList> ParameterStackType;
   typedef symbols<complex>                    ArgumentType;

   static constants                  constants_p;
   static unary_funcs<ElementType>  unary_funcs_p;
   static binary_funcs<ElementType> binary_funcs_p;

   UnitCellParser(ElemStackType& eval_,
                  UnaryFuncStackType& func_stack_,
                  BinaryFuncStackType& bin_func_stack_,
                  IdentifierStackType& IdentifierStack_,
                  NumParameterStackType& NumParameterStack_,
                  FunctionStackType& Functions_,
                  ParameterStackType& Parameters_,
                  ArgumentType& Arguments_,
                  UnitCell const& Cell_, int NumCells_)
      : eval(eval_), func_stack(func_stack_), bin_func_stack(bin_func_stack_),
        IdentifierStack(IdentifierStack_), NumParameterStack(NumParameterStack_),
        FunctionStack(Functions_),
        ParameterStack(Parameters_), Arguments(Arguments_),
        Cell(Cell_), NumCells(NumCells_)
   {}

   template <typename ScannerT>
   struct definition
   {
      definition(UnitCellParser const& self)
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

         filename = lexeme_d[*(anychar_p - chset<>(","))]
            [push_identifier(self.IdentifierStack)];

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

         num_cells_t = (eps_p((str_p("cells") | str_p("sites")) >> '=')
                      >> ((str_p("cells") >> '=' >> expression_t >> ',')
                          | (str_p("sites") >> '=' >> expression_t >> ',')));

         num_cells = (eps_p((str_p("cells") | str_p("sites")) >> '=')
                      >> ((str_p("cells") >> '=' >> expression >> ',')
                          [scale_cells_to_sites(self.Cell, self.NumCells, self.eval)]
                          | (str_p("sites") >> '=' >> expression >> ',')))
            | eps_p[push_number_of_sites(self.Cell, self.NumCells, self.eval)];

         coarse_grain_expression_t = str_p("coarse_grain")
            >> '('
            >> (num_cells_t
                >> expression >> ')');

         coarse_grain_expression = str_p("coarse_grain")
            >> '('
            >> (num_cells
                >> expression >> ')')[push_coarse_grain(self.Cell, self.NumCells, self.eval)];

         fsup_t = str_p("fsup") >> '(' >> expression_t >> ',' >> expression_t >> ',' >> expression_string >> ')';

         fsup = str_p("fsup") >> '(' >> expression >> ',' >> expression >> ','
                              >> expression_string[push_fsup(self.Cell, self.NumCells, self.eval)]
                              >> ')';

         randexpr_t = str_p("rand") >> !('[' >> expression_t >> ',' >> expression_t >> ']')
                                    >> '(' >> expression_t >> *(',' >> expression_t) >> ')';

         randexpr = (str_p("rand") >> (('[' >> expression >> ',' >> expression >> ']')
                                       | eps_p[push_default_randexpr(self.eval)])
                     >> '(' >> expression[init_num_param_stack(self.NumParameterStack)]
                     >> *(',' >> expression[increment_num_param_stack(self.NumParameterStack)]) >> ')')
            [push_rand_expr(self.NumParameterStack, self.eval)];

         swapb_cell_expr_t = (str_p("swapb") >> '(' >> expression_t >> ',' >> expression_t >> ')')
            >> !('[' >> expression_t >> ',' >> expression_t >> ']');

         swapb_cell_expr = (str_p("swapb") >> '(' >> expression >> ',' >> expression >> ')')
            >> (('[' >> expression >> ',' >> expression >> ']')
                [push_swapb_cell_site(self.Cell, self.NumCells, self.eval)]
                |  eps_p[push_swapb_cell(self.Cell, self.NumCells, self.eval)]);

         swapb_site_expr_t = (str_p("swapb") >> '[' >> expression_t >> ',' >> expression_t >> ']');

         swapb_site_expr = (str_p("swapb") >> '[' >> expression >> ',' >> expression >> ']')
            [push_swapb_site(self.Cell, self.NumCells, self.eval)];

         swap_cell_expr_t = (str_p("swap") >> '(' >> expression_t >> ',' >> expression_t >> ')')
            >> !('[' >> expression_t >> ',' >> expression_t >> ']');

         swap_cell_expr = (str_p("swap") >> '(' >> expression >> ',' >> expression >> ')')
            >> (('[' >> expression >> ',' >> expression >> ']')
                [push_swap_cell_site(self.Cell, self.NumCells, self.eval)]
                |  eps_p[push_swap_cell(self.Cell, self.NumCells, self.eval)]);

         swap_site_expr_t = (str_p("swap") >> '[' >> expression_t >> ',' >> expression_t >> ']');

         swap_site_expr =  (str_p("swap") >> '[' >> expression >> ',' >> expression >> ']')
            [push_swap_site(self.Cell, self.NumCells, self.eval)];

         // A local function of the form f(Cell)[Site]{Params}
         local_function_t = identifier_t >> bracket_expr_t >> sq_bracket_expr_t
                                         >> parameter_list_t;

         local_function = identifier[push_function(self.FunctionStack,
                                                   self.ParameterStack)]
            >> bracket_expr
            >> sq_bracket_expr
            >> parameter_list[eval_local_function(self.Cell,
                                                  self.FunctionStack,
                                                  self.ParameterStack,
                                                  self.eval)];

         local_c_function_t = identifier_t >> sq_bracket_expr_t >> parameter_list_t;

         // We can parse a function of the form f[site]{Params}
         // (ie. no cell index), as long as it evalates to a c-number.
         local_c_function = identifier[push_function(self.FunctionStack,
                                                     self.ParameterStack)]
            >> sq_bracket_expr
            >> parameter_list[eval_local_c_function(self.Cell,
                                                    self.FunctionStack,
                                                    self.ParameterStack,
                                                    self.eval)];

         // A local operator of the form Op(cell)[site]
         local_operator_t = identifier_t >> bracket_expr_t >> sq_bracket_expr_t;

         local_operator = identifier[push_identifier(self.IdentifierStack)]
            >> bracket_expr
            >> sq_bracket_expr[push_local_operator(self.Cell,
                                                   self.NumCells,
                                                   self.IdentifierStack,
                                                   self.eval)];

         // A cell function of the form Op(cell){params}
         cell_function_t = identifier_t >> bracket_expr_t >> parameter_list_t;

         cell_function = identifier[push_function(self.FunctionStack, self.ParameterStack)]
            >> bracket_expr
            >> parameter_list[eval_cell_function(self.Cell,
                                                 self.NumCells,
                                                 self.FunctionStack,
                                                 self.ParameterStack,
                                                 self.eval)];

         // A cell operator of the form Op(cell)
         cell_operator_t = identifier_t >> bracket_expr_t;

         cell_operator = identifier[push_identifier(self.IdentifierStack)]
            >> bracket_expr[push_cell_operator(self.Cell,
                                               self.NumCells,
                                               self.IdentifierStack,
                                               self.eval)];

         // An operator of the form Op[i]
         // with no cell index, this can only be a local argument
         local_arg_t = identifier_t >> sq_bracket_expr_t;

         local_arg = identifier[push_identifier(self.IdentifierStack)]
            >> sq_bracket_expr[push_local_argument(self.Cell,
                                                   self.IdentifierStack,
                                                   self.eval)];

         // An operator with no indices, f{Params},
         // is a valid function call, but only if it is c-number valued.
         cell_c_function_t = identifier_t >> parameter_list_t;

         cell_c_function = identifier[push_function(self.FunctionStack,
                                                    self.ParameterStack)]
            >> parameter_list[eval_cell_c_function(self.Cell,
                                                   self.FunctionStack,
                                                   self.ParameterStack,
                                                   self.eval)];

         // combining the above operator/function expressions.
         // It is important that we parse them in the right order, from
         // longest to shortest
         operator_expression_t =
            local_function_t
            | local_c_function_t
            | local_operator_t
            | cell_function_t
            | cell_operator_t
            | local_arg_t
            | cell_c_function_t;

         operator_expression =
            (eps_p(local_function_t) >> local_function)
            | (eps_p(local_c_function_t) >> local_c_function)
            | (eps_p(local_operator_t) >> local_operator)
            | (eps_p(cell_function_t) >> cell_function)
            | (eps_p(cell_operator_t) >> cell_operator)
            | (eps_p(local_arg_t) >> local_arg)
            | (eps_p(cell_c_function_t) >> cell_c_function)
            ;

         expression_string_t = lexeme_d[+((anychar_p - chset<>("()"))
                                          | (ch_p('(') >> expression_string >> ch_p(')')))];

         expression_string = lexeme_d[+((anychar_p - chset<>("()"))
                                        | (ch_p('(') >> expression_string >> ch_p(')')))];

         string_expression_t = str_p("string")
            >> '('
            >> expression_string_t
            >> ')';

         string_expression = str_p("string")
            >> '('
            >> expression_string[push_string(self.Cell, self.NumCells, self.eval)]
            >> ')';

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

         filegrid_expression = str_p("filegrid")
            >> '('
            >> filename >> ','
            >> expression[init_num_param_stack(self.NumParameterStack)]
            >> (*(',' >> expression[increment_num_param_stack(self.NumParameterStack)]) >> ')')
               [eval_filegrid<ElementType, UnitCell>(self.Cell, self.IdentifierStack, self.NumParameterStack, self.eval)];

         filegrid_expression_t = str_p("filegrid")
            >> '('
            >> filename >> ','
            >> expression
            >> *(',' >> expression[increment_num_param_stack(self.NumParameterStack)]) >> ')';

         factor_t =
            imag_t
            |   real_t
            |   unary_function_t
            |   binary_function_t
            |   keyword_d[constants_p[eps_p]]
            |   keyword_d[self.Arguments]
            |   prod_expression_t
            |   commutator_bracket_t
            |   fsup_t
            |   randexpr_t
            |   coarse_grain_expression_t
            |   swapb_cell_expr_t
            |   swapb_site_expr_t
            |   swap_cell_expr_t
            |   swap_site_expr_t
            |   filegrid_expression_t
            |   string_expression_t
            |   '(' >> expression_t >> ')'
            |   ('-' >> factor_t)
            |   ('+' >> factor_t)
            |   operator_expression_t
            ;

         factor =
            imag
            |   real
            |   unary_function
            |   binary_function
            |   keyword_d[constants_p[push_value<ElementType>(self.eval)]]
            |   keyword_d[self.Arguments[push_value<ElementType>(self.eval)]]
            |   prod_expression
            |   commutator_bracket
            |   fsup
            |   randexpr
	    |   coarse_grain_expression
            |   swapb_cell_expr
            |   swapb_site_expr
            |   swap_cell_expr
            |   swap_site_expr
            |   filegrid_expression
            |   string_expression
            |   '(' >> expression >> ')'
            |   ('-' >> factor)[do_negate<ElementType>(self.eval)]
            |   ('+' >> factor)
            |   operator_expression
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
         swap_cell_expr, swap_site_expr, swapb_cell_expr, swapb_site_expr,
         local_function, local_c_function, local_operator, cell_function, cell_operator,
         local_arg, cell_c_function,
         operator_expression, fsup, randexpr,
         expression_string, filename, filegrid_expression,
         string_expression,
         unary_function, binary_function,
         commutator_bracket,
         factor, pow_term, term, expression, num_cells, coarse_grain_expression;

      rule<ScannerT> real_t, imag_t, identifier_t, quantumnumber_t,
         parameter_t, named_parameter_t, parameter_list_t,
         prod_expression_t,
         bracket_expr_t, sq_bracket_expr_t,
         swap_cell_expr_t, swap_site_expr_t, swapb_cell_expr_t, swapb_site_expr_t,
         local_function_t, local_c_function_t, local_operator_t, cell_function_t, cell_operator_t,
         local_arg_t, cell_c_function_t,
         operator_expression_t, fsup_t, randexpr_t,
         expression_string_t,
         string_expression_t, filegrid_expression_t,
         unary_function_t, binary_function_t,
         commutator_bracket_t,
         factor_t, pow_term_t, term_t, expression_t, num_cells_t, coarse_grain_expression_t;

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
   UnitCell const& Cell;
   int NumCells;
};


// global variables (static members of UnitCellParser)
constants UnitCellParser::constants_p;
unary_funcs<UnitCellParser::ElementType> UnitCellParser::unary_funcs_p;
binary_funcs<UnitCellParser::ElementType> UnitCellParser::binary_funcs_p;

UnitCellElementType
ParseUnitCellElement(UnitCell const& Cell, int NumCells, std::string const& Str,
                     Function::ArgumentList const& Args)
{
   typedef UnitCellParser::ElementType ElementType;

   UnitCellParser::ElemStackType         ElemStack;
   UnitCellParser::UnaryFuncStackType    UnaryFuncStack;
   UnitCellParser::BinaryFuncStackType   BinaryFuncStack;
   UnitCellParser::IdentifierStackType   IdentStack;
   UnitCellParser::NumParameterStackType NumParamStack;
   UnitCellParser::ParameterStackType    ParamStack;
   UnitCellParser::FunctionStackType     FunctionStack;
   UnitCellParser::ArgumentType          Arguments;

   for (Function::ArgumentList::const_iterator I = Args.begin(); I != Args.end(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   // Put the lattice args into Arguments - this won't override existing values
   for (UnitCell::const_argument_iterator I = Cell.begin_arg(); I != Cell.end_arg(); ++I)
   {
      Arguments.add(I->first.c_str(), I->second);
   }

   char const* beg = Str.c_str();
   char const* end = beg + Str.size();

   UnitCellParser Parser(ElemStack, UnaryFuncStack, BinaryFuncStack, IdentStack, NumParamStack,
                         FunctionStack, ParamStack,
                         Arguments, Cell, NumCells);

   try
   {
      CheckParentheses(beg, end);

      parse_info<> info = parse(beg, Parser, space_p);
      if (!info.full)
         throw ParserError::AtRange("Failed to parse an expression", info.stop, end);
   }
   catch (ParserError const& p)
   {
      throw ParserError::Finalize(p, "while parsing a unit cell operator:", beg, end);
   }
   catch (std::exception const& p)
   {
      throw ParserError::Finalize(p, "while parsing a unit cell operator:", beg, end);
   }
   catch (...)
   {
      throw;
   }

   CHECK(UnaryFuncStack.empty());
   CHECK(BinaryFuncStack.empty());
   CHECK(IdentStack.empty());
   CHECK(ParamStack.empty());
   CHECK(NumParamStack.empty());
   CHECK(!ElemStack.empty());
   ElementType Result = ElemStack.top();
   ElemStack.pop();
   CHECK(ElemStack.empty())(ElemStack.size());

   return Result;
}

UnitCellMPO
ParseUnitCellOperator(UnitCell const& Cell, int NumCells, std::string const& Str,
                      Function::ArgumentList const& Args, InfiniteLattice const* Lat)
{
   InfiniteLattice const* OldILattice = ILattice;
   if (Lat)
   {
      ILattice = Lat;
   }
   ElementType Result = ParseUnitCellElement(Cell, NumCells, Str, Args);
   if (Lat)
   {
      ILattice = OldILattice;
   }

   UnitCellMPO* Op = boost::get<UnitCellMPO>(&Result);
   if (Op)
   {
      return *Op;
   }
   // else, we also handle the case where the operator is a number
   complex x = boost::get<complex>(Result);

   return x*Cell["I"];
}

UnitCellElementType
ParseOperator(UnitCell const& Cell, std::string const& Str,
              Function::ArgumentList const& Args)
{
   return ParseUnitCellElement(Cell, 0, Str, Args);
}

std::complex<double>
ParseUnitCellNumber(UnitCell const& Cell, int NumCells, std::string const& Str,
                    Function::ArgumentList const& Args)
{
   UnitCellElementType Result = ParseUnitCellElement(Cell, NumCells, Str, Args);
   complex* x = boost::get<complex>(&Result);
   CHECK(x)("expression is an operator, but a c-number was expected!")(Str);
   return *x;
}

std::pair<UnitCellMPO, InfiniteLattice>
ParseUnitCellOperatorAndLattice(std::string const& Str)
{
   std::string::const_iterator Delim = std::find(Str.begin(), Str.end(), ':');
   if (Delim == Str.end())
   {
      PANIC("fatal: expression of the form \"lattice:expression\" expected.")(Str);
   }

   std::string LatticeFile = std::string(Str.begin(), Delim);
   boost::trim(LatticeFile);
   pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

   InfiniteLattice const* OldILattice = ILattice;
   ILattice = &(*Lattice);

   ++Delim;
   std::string Expr(Delim, Str.end());

   UnitCellMPO Op = ParseUnitCellOperator(Lattice->GetUnitCell(), 0, Expr, Lattice->args());

   ILattice = OldILattice;

   return std::make_pair(Op, *Lattice);
}
