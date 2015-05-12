// -*- C++ -*- $Id$

#include "unitcell-parser.h"
#include "parser/parser.h"
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>

using namespace Parser;

typedef boost::variant<complex, UnitCellMPO> element_type;

namespace UP
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

// We have 5 variants of operator expressions.
// Op(n)[j]  - local operator at site j of the n'th unit cell
// Op[j](n)  - local operator at site j of the n'th unit cell
// Op[j]     - requires that there is exactly 1 unit cell, the unit cell number is implicitly 0
// Op(n)     - a unit cell operator on the n'th unit cell
// Op        - requires that there is exactly 1 unit cell, the unit cell number is implicitly 0
//
// The same 5 variants exist also with a parameter pack

// a unit cell operator that has no unit cell index attached.
// This is only possible if we only have one cell (NumCells == 1)
struct push_operator
{
   push_operator(UnitCell const& Cell_, int NumCells_,
		 std::stack<std::string>& IdentStack_, 
		 std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      CHECK_EQUAL(NumCells, 1)("Operator must supply a cell index")(OpName);
      CHECK(Cell.operator_exists(OpName))("Operator does not exist in the unit cell")(OpName);

      eval.push(element_type(Cell.Operator(OpName)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<element_type >& eval;
};

// a unit cell operator with cell index attached.
// The cell number j must satisfy j >= 0 && (NumCells == 0 || j < NumCells)
struct push_operator_cell
{
   push_operator_cell(UnitCell const& Cell_, int NumCells_,
		      std::stack<std::string>& IdentStack_, 
		      std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Cell index out of bounds")(j)(NumCells);
      CHECK(Cell.operator_exists(OpName))("Operator does not exist in the unit cell")(OpName);

      eval.push(element_type(Cell.OperatorAtCell(OpName, j)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<element_type >& eval;
};

// a local operator that has no cell number, only a site attached
// This is only possible if we only have one cell (NumCells == 1)
struct push_operator_site
{
   push_operator_site(UnitCell const& Cell_, int NumCells_,
		      std::stack<std::string>& IdentStack_, 
		      std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      CHECK_EQUAL(NumCells, 1)("Operator must supply a cell index")(OpName);

      int n = pop_int(eval);
      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);
      eval.push(element_type(Cell.LocalOperator(OpName, n)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<element_type >& eval;
};

// a local operator that specifies a cell index as well.
// The cell index was pushed first, in the notation Op(n)[j]
struct push_operator_cell_site
{
   push_operator_cell_site(UnitCell const& Cell_, int NumCells_,
			   std::stack<std::string>& IdentStack_, 
			   std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int n = pop_int(eval);
      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Site index out of bounds")(j)(NumCells);
      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);

      // Fetch the operator and JW string
      eval.push(element_type(Cell.LocalOperator(OpName, j, n)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<element_type >& eval;
};

// a local operator that specifies a cell index as well.
// The site number was pushed first, in the notation Op[j](n)
struct push_operator_site_cell
{
   push_operator_site_cell(UnitCell const& Cell_, int NumCells_,
			   std::stack<std::string>& IdentStack_, 
			   std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Cell index out of bounds")(j)(NumCells);

      int n = pop_int(eval);
      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);

      // Fetch the operator and JW string
      eval.push(element_type(Cell.LocalOperator(OpName, j, n)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<element_type >& eval;
};

//
// operators with parameters
//

struct push_operator_param
{
   push_operator_param(UnitCell const& Cell_, int NumCells_,
		       std::stack<std::string>& IdentStack_, 
		       std::stack<std::stack<element_type> >& ParamStack_,
		       std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), 
	ParamStack(ParamStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      CHECK_EQUAL(NumCells, 1)("Operator must supply a cell index")(OpName);

      // Convert the parameter stack into a vector of numbers
      int NumParam = ParamStack.top().size(); // number of parameters
      std::vector<std::complex<double> > Params(NumParam);
      while (!ParamStack.top().empty())
      {
	 Params[--NumParam] = boost::get<std::complex<double> >(ParamStack.top().top());
	 ParamStack.top().pop();
      }
      ParamStack.pop();

      eval.push(element_type(Cell.OperatorFunction(OpName, Params)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<std::stack<element_type> >& ParamStack;
   std::stack<element_type >& eval;
};

// a unit cell operator with cell index attached.
// The cell number j must satisfy j >= 0 && (NumCells == 0 || j < NumCells)
struct push_operator_cell_param
{
   push_operator_cell_param(UnitCell const& Cell_, int NumCells_,
			    std::stack<std::string>& IdentStack_, 
			    std::stack<std::stack<element_type> >& ParamStack_,
			    std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_),
	ParamStack(ParamStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Cell index out of bounds")(j)(NumCells);
 
      // Convert the parameter stack into a vector of numbers
      int NumParam = ParamStack.top().size(); // number of parameters
      std::vector<std::complex<double> > Params(NumParam);
      while (!ParamStack.top().empty())
      {
	 Params[--NumParam] = boost::get<std::complex<double> >(ParamStack.top().top());
	 ParamStack.top().pop();
      }
      ParamStack.pop();

      // Construct the operator
      eval.push(element_type(Cell.OperatorFunction(OpName, Params)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<std::stack<element_type> >& ParamStack;
   std::stack<element_type >& eval;
};

// a local operator that has no cell number, only a site attached
// This is only possible if we only have one cell (NumCells == 1)
struct push_operator_site_param
{
   push_operator_site_param(UnitCell const& Cell_, int NumCells_,
			    std::stack<std::string>& IdentStack_, 
			    std::stack<std::stack<element_type> >& ParamStack_,
			    std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_),
	ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      CHECK_EQUAL(NumCells, 1)("Operator must supply a cell index")(OpName);

      int n = pop_int(eval);

      // Convert the parameter stack into a vector of numbers
      int NumParam = ParamStack.top().size(); // number of parameters
      std::vector<std::complex<double> > Params(NumParam);
      while (!ParamStack.top().empty())
      {
	 Params[--NumParam] = boost::get<std::complex<double> >(ParamStack.top().top());
	 ParamStack.top().pop();
      }
      ParamStack.pop();

      eval.push(element_type(Cell.OperatorFunction(OpName, n, Params)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<std::stack<element_type> >& ParamStack;
   std::stack<element_type >& eval;
};

// a local operator that specifies a cell index as well.
// The cell index was pushed first, in the notation Op(n)[j]
struct push_operator_cell_site_param
{
   push_operator_cell_site_param(UnitCell const& Cell_, int NumCells_,
				 std::stack<std::string>& IdentStack_, 
				 std::stack<std::stack<element_type> >& ParamStack_,
				 std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), 
	ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int n = pop_int(eval);

      complex CellIndex = boost::get<complex>(eval.top());
      eval.pop();
      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Cell index out of bounds")(j)(NumCells);

      // Convert the parameter stack into a vector of numbers
      int NumParam = ParamStack.top().size(); // number of parameters
      std::vector<std::complex<double> > Params(NumParam);
      while (!ParamStack.top().empty())
      {
	 Params[--NumParam] = boost::get<std::complex<double> >(ParamStack.top().top());
	 ParamStack.top().pop();
      }
      ParamStack.pop();
      eval.push(element_type(Cell.OperatorFunction(OpName, n, Params)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentStack;
   std::stack<std::stack<element_type> >& ParamStack;
   std::stack<element_type >& eval;
};

// a local operator that specifies a cell index as well.
// The site number was pushed first, in the notation Op[j](n)
struct push_operator_site_cell_param
{
   push_operator_site_cell_param(UnitCell const& Cell_, int NumCells_,
				 std::stack<std::string>& IdentStack_, 
				 std::stack<std::stack<element_type> >& ParamStack_,
				 std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentStack(IdentStack_), 
	ParamStack(ParamStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentStack.top();
      IdentStack.pop();

      int j = pop_int(eval);
      CHECK(NumCells == 0 || (j >= 0 && j < NumCells))("Cell index out of bounds")(j)(NumCells);

      int n = pop_int(eval);

      // Convert the parameter stack into a vector of numbers
      int NumParam = ParamStack.top().size(); // number of parameters
      std::vector<std::complex<double> > Params(NumParam);
      while (!ParamStack.top().empty())
      {
	 Params[--NumParam] = boost::get<std::complex<double> >(ParamStack.top().top());
	 ParamStack.top().pop();
      }
      ParamStack.pop();

      // Fetch the operator and JW string
      eval.push(element_type(Cell.OperatorFunction(OpName, n, Params)));
   }

   UnitCell const& Cell;
   int NumCells;
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

//
// swapb - swap without fermion signs
//

// notation for swapb(cell1,cell2)[site1,site2]
struct push_swapb_cell_site
{
   push_swapb_cell_site(UnitCell const& Cell_, int NumCells_,
		       std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);
      int Cell2 = pop_int(eval);
      int Cell1 = pop_int(eval);

      CHECK(NumCells == 0 || (Cell1 >= 0 && Cell1 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);
      CHECK(NumCells == 0 || (Cell2 >= 0 && Cell2 < NumCells))("Cell index out of bounds")(Cell1)(NumCells);

      eval.push(element_type(Cell.swap_gate_no_sign(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

// notation for swapb[site1,site2]
struct push_swapb_site
{
   push_swapb_site(UnitCell const& Cell_, int NumCells_,
		  std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      CHECK(NumCells == 0 || NumCells == 1)("Cell index required");
      int Cell2 = 0;
      int Cell1 = 0;
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);

      eval.push(element_type(Cell.swap_gate_no_sign(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

// notation for swapb(cell1,cell2)
struct push_swapb_cell
{
   push_swapb_cell(UnitCell const& Cell_, int NumCells_,
		  std::stack<element_type>& eval_)
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

      eval.push(element_type(Cell.swap_gate_no_sign(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

//
// swap (with fermion signs)
//

// notation for swap(cell1,cell2)[site1,site2]
struct push_swap_cell_site
{
   push_swap_cell_site(UnitCell const& Cell_, int NumCells_,
		       std::stack<element_type>& eval_)
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
      //eval.push(element_type(Cell.swap_gate(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

// notation for swap[site1,site2]
struct push_swap_site
{
   push_swap_site(UnitCell const& Cell_, int NumCells_,
		  std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      CHECK(NumCells == 0 || NumCells == 1)("Cell index required");
      int Cell2 = 0;
      int Cell1 = 0;
      int Site2 = pop_int(eval);
      int Site1 = pop_int(eval);

      PANIC("Swap() is not yet implemented");
      //eval.push(element_type(Cell.swap_gate(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

// notation for swapb(cell1,cell2)
struct push_swap_cell
{
   push_swap_cell(UnitCell const& Cell_, int NumCells_,
		  std::stack<element_type>& eval_)
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
      //      eval.push(element_type(Cell.swap_gate(Cell1, Site1, Cell2, Site2)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

struct push_string
{
   push_string(UnitCell const& Cell_, int NumCells_,
	       std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), eval(eval_) 
   {
      // If the number of unit cells is unspecified, use 1
      if (NumCells == 0)
	 NumCells = 1;
   }
   
   void operator()(char const* Start, char const* End) const
   {
      FiniteMPO Op = ParseStringOperator(*Cell.GetSiteList(), 
					 std::string(Start, End), NumCells*Cell.size());
      eval.push(UnitCellMPO(Cell.GetSiteList(), Op, LatticeCommute::Bosonic, 0));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<element_type >& eval;
};

} // namespace UP

using namespace UP;

struct UnitCellParser : public grammar<UnitCellParser>
{
   typedef boost::variant<complex, UnitCellMPO> element_type;
   typedef boost::function<element_type(element_type)> unary_func_type;
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;

   typedef std::stack<element_type>       ElemStackType;
   typedef std::stack<std::string>        IdentStackType;
   typedef std::stack<unary_func_type>    UnaryFuncStackType;
   typedef std::stack<binary_func_type>   BinaryFuncStackType;
   typedef std::stack<ElemStackType>      ParamStackType;

   static constants                  constants_p;
   static unary_funcs<element_type>  unary_funcs_p;
   static binary_funcs<element_type> binary_funcs_p;

   UnitCellParser(ElemStackType& eval_, 
		  IdentStackType& identifier_stack_,
		  ParamStackType& param_stack_,
		  UnaryFuncStackType& func_stack_,
		  BinaryFuncStackType& bin_func_stack_,
		  UnitCell const& Cell_, int NumCells_)
      : eval(eval_), identifier_stack(identifier_stack_), 
	param_stack(param_stack_),
	func_stack(func_stack_), bin_func_stack(bin_func_stack_), Cell(Cell_), NumCells(NumCells_)
   {}
   
   template <typename ScannerT>
   struct definition
   {
      definition(UnitCellParser const& self)
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

	 swapb_cell_expr = (str_p("swapb") >> '(' >> expression >> ',' >> expression >> ')')
	    >> (('[' >> expression >> ',' >> expression >> ']')[push_swapb_cell_site(self.Cell, self.NumCells, self.eval)]
		|  eps_p[push_swapb_cell(self.Cell, self.NumCells, self.eval)]);

	 swapb_site_expr =  (str_p("swapb") >> '[' >> expression >> ',' >> expression >> ']')
	    [push_swapb_site(self.Cell, self.NumCells, self.eval)];

	 swap_cell_expr = (str_p("swap") >> '(' >> expression >> ',' >> expression >> ')')
	    >> (('[' >> expression >> ',' >> expression >> ']')[push_swap_cell_site(self.Cell, self.NumCells, self.eval)]
		|  eps_p[push_swap_cell(self.Cell, self.NumCells, self.eval)]);

	 swap_site_expr =  (str_p("swap") >> '[' >> expression >> ',' >> expression >> ']')
	    [push_swap_site(self.Cell, self.NumCells, self.eval)];

	 // an operator of the form Op(u)[i]
	 // If there is a parameter pack then it is a local operator function
	 // otherwise it is a local operator
	 operator_bracket_sq = sq_bracket_expr
	    >> (parameter_list[push_operator_cell_site_param(self.Cell, 
							     self.NumCells, 
							     self.identifier_stack, 
							     self.param_stack,
							     self.eval)]
		|   eps_p[push_operator_cell_site(self.Cell, 
						  self.NumCells, 
						  self.identifier_stack, 
						  self.eval)]
		);

	 // An operator of the form Op(u)
	 // This could be followed by:
	 // square brackets - a local operator
	 // a parameter list - a unit cell operator function
	 // otherwise it is a unit cell operator
	 operator_bracket = bracket_expr
	    >> (operator_bracket_sq
		|   parameter_list[push_operator_cell_param(self.Cell, 
							    self.NumCells, 
							    self.identifier_stack, 
							    self.param_stack,
							    self.eval)]
		|   eps_p[push_operator_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
		);

	 // An operator of the form Op[i](a)
	 // If there is a parameter pack then it is a local operator function
	 // otherwise it is a local operator
	 operator_sq_bracket = bracket_expr
	    >> (parameter_list[push_operator_site_cell_param(self.Cell, 
							    self.NumCells, 
							    self.identifier_stack, 
							    self.param_stack,
							    self.eval)]
		|   eps_p[push_operator_site_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
		);

	 // An operator of the form Op[i]
	 // This could be followed by:
	 // brackets - a local operator
	 // a parameter list - a unit cell operator function with implicit unit cell number
	 // otherwise it is a local operator with implicit unit cell number
	 operator_sq = sq_bracket_expr
	    >> (operator_sq_bracket
		|   parameter_list[push_operator_site_param(self.Cell, 
							    self.NumCells, 
							    self.identifier_stack, 
							    self.param_stack,
							    self.eval)]
		| eps_p[push_operator_site(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
		);

	 operator_expression = 
	    identifier
	    >> (operator_bracket
		|   operator_sq_bracket
		|   parameter_list[push_operator_param(self.Cell, 
						       self.NumCells, 
						       self.identifier_stack, 
						       self.param_stack,
						       self.eval)]
		|   eps_p
       [push_operator(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
		);
	 
	 expression_string = lexeme_d[+((anychar_p - chset<>("()"))
					| (ch_p('(') >> expression_string >> ch_p(')')))];

	 string_expression = str_p("string")
	    >> '(' 
	    >> expression_string[push_string(self.Cell, self.NumCells, self.eval)]
	    >> ')';

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
	    |   commutator_bracket
	    |   swapb_cell_expr
	    |   swapb_site_expr
	    |   swap_cell_expr
	    |   swap_site_expr
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
	 binary_function, bracket_expr, quantumnumber, prod_expression, sq_bracket_expr, 
	 operator_expression, operator_bracket_sq, operator_sq_bracket, operator_bracket, operator_sq,
	 parameter, parameter_list, swap_cell_expr, swap_site_expr, swapb_cell_expr, swapb_site_expr,
	 local_operator, local_operator_cell_site, local_operator_site_cell, 
	 expression_string, string_expression,
	 identifier, pow_term, commutator_bracket;
      rule<ScannerT> const&
      start() const { return expression; }
   };
   
   std::stack<element_type>& eval;
   IdentStackType& identifier_stack;
   ParamStackType& param_stack;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   UnitCell const& Cell;
   int NumCells;
};


// global variables (static members of UnitCellParser)
constants UnitCellParser::constants_p;
unary_funcs<UnitCellParser::element_type> UnitCellParser::unary_funcs_p;
binary_funcs<UnitCellParser::element_type> UnitCellParser::binary_funcs_p;

UnitCellMPO
ParseUnitCellOperator(UnitCell const& Cell, int NumCells, std::string const& Str)
{
   typedef UnitCellParser::element_type element_type;

   UnitCellParser::ElemStackType ElemStack;
   UnitCellParser::IdentStackType IdentStack;
   UnitCellParser::ParamStackType ParamStack;
   UnitCellParser::UnaryFuncStackType UnaryFuncStack;
   UnitCellParser::BinaryFuncStackType BinaryFuncStack;

   UnitCellParser Parser(ElemStack, IdentStack, ParamStack, 
			 UnaryFuncStack, BinaryFuncStack, Cell, NumCells);

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

   UnitCellMPO* Op = boost::get<UnitCellMPO>(&Result);
   if (Op)
   {
      return *Op;
   }
   // else, we also handle the case where the operator is a number
   complex x = boost::get<complex>(Result);
   return x*Cell.Operator("I");
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

   ++Delim;
   std::string Expr(Delim, Str.end());

   UnitCellMPO Op = ParseUnitCellOperator(Lattice->GetUnitCell(), 0, Expr);
   return std::make_pair(Op, *Lattice);
}
