// -*- C++ -*- $Id$

#include "unitcell-parser.h"
#include "parser/parser.h"


using namespace Parser;

typedef boost::variant<complex, FiniteMPO> element_type;

namespace Parser
{

// We have a problem: for operators with indefinite support, we need to extend MPOs as necessary
// so that binary operations act on operands with the same number of unit cells.  We can't do that
// without knowing the actual unit cell, since it might be some operator that mixes the basis around etc.
// This means that we need to specialize the functions from parser.h

template <>
struct binary_addition<element_type> : boost::static_visitor<element_type>
{
   binary_addition(UnitCell const& Cell_) 
      : Cell(Cell_) {}

   element_type operator()(complex const& x, complex const& y) const
   { 
      return element_type(x+y);
   }

   element_type operator()(complex const& x, FiniteMPO const& y) const
   {
      return element_type(x*MakeIdentityFrom(y) + y);
   }

   element_type operator()(FiniteMPO const& x, complex const& y) const
   {
      return element_type(x + y*MakeIdentityFrom(x));
   }

   element_type operator()(FiniteMPO const& x, FiniteMPO const& y) const
   {
      if (x.size() < y.size())
      {
	 FiniteMPO xExtend = join(x, repeat(identity_mpo(Cell, x.qn2()), (y.size()-x.size())/Cell.size()));
	 return element_type(xExtend+y);
      }
      else if (x.size() > y.size())
      {
	 FiniteMPO yExtend = join(y, repeat(identity_mpo(Cell, y.qn2()), (x.size()-y.size())/Cell.size()));
	 return element_type(x+yExtend);
      }
      // else
      return element_type(x+y);
   }

   UnitCell const& Cell;
};

template <>
struct binary_subtraction<element_type> : boost::static_visitor<element_type>
{
   binary_subtraction(UnitCell const& Cell_) 
      : Cell(Cell_) {}

   element_type operator()(complex const& x, complex const& y) const
   { 
      return element_type(x-y);
   }

   element_type operator()(complex const& x, FiniteMPO const& y) const
   {
      return element_type(x*MakeIdentityFrom(y) - y);
   }

   element_type operator()(FiniteMPO const& x, complex const& y) const
   {
      return element_type(x - y*MakeIdentityFrom(x));
   }

   element_type operator()(FiniteMPO const& x, FiniteMPO const& y) const
   {
      if (x.size() < y.size())
      {
	 FiniteMPO xExtend = join(x, repeat(identity_mpo(Cell, x.qn2()), (y.size()-x.size())/Cell.size()));
	 return element_type(xExtend-y);
      }
      else if (x.size() > y.size())
      {
	 FiniteMPO yExtend = join(y, repeat(identity_mpo(Cell, y.qn2()), (x.size()-y.size())/Cell.size()));
	 return element_type(x-yExtend);
      }
      // else
      return element_type(x-y);
   }

   UnitCell const& Cell;
};

template <>
struct binary_commutator<element_type> : boost::static_visitor<element_type>
{
   binary_commutator(UnitCell const& Cell_) 
      : Cell(Cell_) {}

   element_type operator()(complex const& x, complex const& y) const
   { 
      return complex(0.0);
   }

   element_type operator()(complex const& x, FiniteMPO const& y) const
   {
      return complex(0.0);
   }

   element_type operator()(FiniteMPO const& x, complex const& y) const
   {
      return complex(0.0);
   }

   element_type operator()(FiniteMPO const& x, FiniteMPO const& y) const
   {
      if (x.size() < y.size())
      {
	 FiniteMPO xExtend = join(x, repeat(identity_mpo(Cell, x.qn2()), (y.size()-x.size())/Cell.size()));
	 return element_type(xExtend*y - y*xExtend);
      }
      else if (x.size() > y.size())
      {
	 FiniteMPO yExtend = join(y, repeat(identity_mpo(Cell, y.qn2()), (x.size()-y.size())/Cell.size()));
	 return element_type(x*yExtend - yExtend*x);
      }
      // else
      return element_type(x*y - y*x);
   }

   UnitCell const& Cell;
};

template <>
struct binary_multiplication<element_type> : boost::static_visitor<element_type>
{
   binary_multiplication(UnitCell const& Cell_) 
      : Cell(Cell_) {}

   template <typename T1, typename T2>
   element_type operator()(T1 const& x, T2 const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(FiniteMPO const& x, FiniteMPO const& y) const
   {
      if (x.size() < y.size())
      {
	 FiniteMPO xExtend = join(x, repeat(identity_mpo(Cell, x.qn2()), (y.size()-x.size())/Cell.size()));
	 return element_type(xExtend*y);
      }
      else if (x.size() > y.size())
      {
	 FiniteMPO yExtend = join(y, repeat(identity_mpo(Cell, y.qn2()), (x.size()-y.size())/Cell.size()));
	 return element_type(x*yExtend);
      }
      // else
      return element_type(x*y);
   }

   UnitCell const& Cell;
};

template <>
struct binary_dot_product<element_type> : boost::static_visitor<element_type>
{
   binary_dot_product(UnitCell const& Cell_) 
      : Cell(Cell_) {}

   element_type operator()(complex const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(complex const& x, FiniteMPO const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(FiniteMPO const& x, complex const& y) const
   {
      return element_type(x*y);
   }

   element_type operator()(FiniteMPO const& x, FiniteMPO const& y) const
   {
      if (x.size() < y.size())
      {
	 FiniteMPO xExtend = join(x, repeat(identity_mpo(Cell, x.qn2()), (y.size()-x.size())/Cell.size()));
	 return dot(xExtend, y);
      }
      else if (x.size() > y.size())
      {
	 FiniteMPO yExtend = join(y, repeat(identity_mpo(Cell, y.qn2()), (x.size()-y.size())/Cell.size()));
	 return dot(x, yExtend);
      }
      // else
      return dot(x,y);
   }

   UnitCell const& Cell;

};

template <>
struct ternary_product<element_type> : boost::static_visitor<element_type>
{
   ternary_product(UnitCell const& Cell_, std::string const& q_)
      : Cell(Cell_), q(Cell_.GetSymmetryList(), q_) {}

   element_type operator()(complex const& x, complex const& y) const
   {
      CHECK(is_scalar(q));
      return element_type(x*y);
   }

   element_type operator()(complex const& x, FiniteMPO const& y) const
   {
      CHECK_EQUAL(q, y.TransformsAs());
      return element_type(x*y);
   }

   element_type operator()(FiniteMPO const& x, complex const& y) const
   {
      CHECK_EQUAL(q, x.TransformsAs());
      return element_type(x*y);
   }

   element_type operator()(FiniteMPO const& x, FiniteMPO const& y) const
   {
      if (x.size() < y.size())
      {
	 FiniteMPO xExtend = join(x, repeat(identity_mpo(Cell, x.qn2()), (y.size()-x.size())/Cell.size()));
	 return prod(xExtend, y, q);
      }
      else if (x.size() > y.size())
      {
	 FiniteMPO yExtend = join(y, repeat(identity_mpo(Cell, y.qn2()), (x.size()-y.size())/Cell.size()));
	 return prod(x, yExtend, q);
      }
      // else
      return prod(x,y,q);
   }

   UnitCell const& Cell;
   QuantumNumbers::QuantumNumber q;
};

template <>
struct push_prod<element_type>
{
   push_prod(UnitCell const& Cell_, std::stack<std::string>& identifier_stack_, std::stack<element_type >& eval_)
      : Cell(Cell_), identifier_stack(identifier_stack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      element_type op2 = eval.top();
      eval.pop();
      element_type op1 = eval.top();
      eval.pop();
      std::string q = identifier_stack.top();
      identifier_stack.pop();

      eval.push(boost::apply_visitor(ternary_product<element_type>(Cell, q), op1, op2));
   }

   UnitCell const& Cell;
   std::stack<std::string>& identifier_stack;
   std::stack<element_type>& eval;
};

template <>
struct binary_funcs<element_type> : symbols<boost::function<element_type(element_type, element_type)> >
{
   typedef boost::function<element_type(element_type, element_type)> binary_func_type;
   binary_funcs(UnitCell const& Cell)
   {
      this->add.operator()
         ("dot", make_apply_binary_math<element_type>(binary_dot_product<element_type>(Cell)))
	 //         ("outer", make_apply_binary_math<element_type>(binary_outer_product<element_type>()))
	 //	 ("cross", make_apply_binary_math<element_type>(binary_cross_product<element_type>()))
         ;
   }
};

} // namespace Parser


// a unit cell operator that has no unit cell index attached.
// This is only possible if we only have one cell (NumCells == 1)
struct push_operator_no_cell
{
   push_operator_no_cell(UnitCell const& Cell_, int NumCells_,
			 std::stack<std::string>& IdentifierStack_, 
			 std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentifierStack(IdentifierStack_), eval(eval_) {}
   
   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      CHECK_EQUAL(NumCells, 1)("Operator must supply a cell index")(OpName);
      CHECK(Cell.operator_exists(OpName))("Operator does not exist in the unit cell")(OpName);

      eval.push(element_type(Cell.Operator(OpName)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
};

// a unit cell operator with cell index attached.
// The cell number j must satisfy j >= 0 && (NumCells == 0 || j < NumCells)
struct push_operator_cell
{
   push_operator_cell(UnitCell const& Cell_, int NumCells_,
		      std::stack<std::string>& IdentifierStack_, 
		      std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentifierStack(IdentifierStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      complex CellIndex = boost::get<complex>(eval.top());
      eval.pop();
      int j = int(CellIndex.real() + 0.5);
      CHECK(norm_frob(CellIndex - double(j)) < 1E-7)("Cell index must be an integer")(CellIndex);
      CHECK(j >= 0)("Cell index must be non-negative")(j)(OpName);
      CHECK(NumCells == 0 || j < NumCells)("Cell index out of bounds")(j)(NumCells);
      CHECK(Cell.operator_exists(OpName))("Operator does not exist in the unit cell")(OpName);

      FiniteMPO CellOperator = Cell.Operator(OpName);
      FiniteMPO JWOperator = jw_string_mpo(Cell, CellOperator);

      // The actual operator is j copies of the JW string, followed by the cell operator
      FiniteMPO Op = join(repeat(JWOperator, j), CellOperator);

      // extend the operator to have support over NumCells
      if (NumCells != 0 && j < NumCells-1)
      {
	 Op = join(Op, repeat(identity_mpo(Cell), NumCells-j-1));
      }

      eval.push(element_type(Op));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
};

// a local operator that has no cell index attached
// This is only possible if we only have one cell (NumCells == 1)
struct push_local_operator_no_cell
{
   push_local_operator_no_cell(UnitCell const& Cell_, int NumCells_,
			       std::stack<std::string>& IdentifierStack_, 
			       std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentifierStack(IdentifierStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      CHECK_EQUAL(NumCells, 1)("Operator must supply a cell index")(OpName);

      complex SiteIndex = boost::get<complex>(eval.top());
      eval.pop();
      int n = int(SiteIndex.real() + 0.5);
      CHECK(norm_frob(SiteIndex - double(n)) < 1E-7)("Site index must be an integer")(SiteIndex);

      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);
      eval.push(element_type(Cell.Operator(OpName, n)));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
};

// a local operator that specifies a cell index as well.
// The cell index was pushed first, in the notation Op(n)[j]
struct push_local_operator_cell_index
{
   push_local_operator_cell_index(UnitCell const& Cell_, int NumCells_,
				  std::stack<std::string>& IdentifierStack_, 
				  std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentifierStack(IdentifierStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      complex SiteIndex = boost::get<complex>(eval.top());
      eval.pop();
      int n = int(SiteIndex.real() + 0.5);
      CHECK(norm_frob(SiteIndex - double(n)) < 1E-7)("Site index must be an integer")(SiteIndex);

      complex CellIndex = boost::get<complex>(eval.top());
      eval.pop();
      int j = int(CellIndex.real() + 0.5);
      CHECK(norm_frob(CellIndex - double(j)) < 1E-7)("Cell index must be an integer")(CellIndex);
      CHECK(j >= 0)("Cell index must be non-negative")(j)(OpName);
      CHECK(NumCells == 0 || j < NumCells)("Cell index out of bounds")(j)(NumCells);
      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);

      // Fetch the operator and JW string
      FiniteMPO CellOperator = Cell.Operator(OpName, n);
      FiniteMPO JWOperator = jw_string_mpo(Cell, CellOperator);

      // The actual operator is j copies of the JW string, followed by the cell operator
      FiniteMPO Op = join(repeat(JWOperator, j), CellOperator);

      // extend the operator to have support over NumCells
      if (NumCells != 0 && j < NumCells-1)
      {
	 Op = join(Op, repeat(identity_mpo(Cell), NumCells-j-1));
      }
      eval.push(element_type(Op));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
};

// a local operator that specifies a cell index as well.
// The site number was pushed first, in the notation Op[j](n)
struct push_local_operator_index_cell
{
   push_local_operator_index_cell(UnitCell const& Cell_, int NumCells_,
				  std::stack<std::string>& IdentifierStack_, 
				  std::stack<element_type>& eval_)
      : Cell(Cell_), NumCells(NumCells_), IdentifierStack(IdentifierStack_), eval(eval_) {}

   void operator()(char const*, char const*) const
   {
      std::string OpName = IdentifierStack.top();
      IdentifierStack.pop();

      complex CellIndex = boost::get<complex>(eval.top());
      eval.pop();
      int j = int(CellIndex.real() + 0.5);
      CHECK(norm_frob(CellIndex - double(j)) < 1E-7)("Cell index must be an integer")(CellIndex);
      CHECK(j >= 0)("Cell index must be non-negative")(j)(OpName);
      CHECK(NumCells == 0 || j < NumCells)("Cell index out of bounds")(j)(NumCells);

      complex SiteIndex = boost::get<complex>(eval.top());
      eval.pop();
      int n = int(SiteIndex.real() + 0.5);
      CHECK(norm_frob(SiteIndex - double(n)) < 1E-7)("Site index must be an integer")(SiteIndex);
      CHECK(Cell.operator_exists(OpName, n))("Local operator does not exist in the unit cell")(OpName)(n);

      // Fetch the operator and JW string
      FiniteMPO CellOperator = Cell.Operator(OpName, n);
      FiniteMPO JWOperator = jw_string_mpo(Cell, CellOperator);

      // The actual operator is j copies of the JW string, followed by the cell operator
      FiniteMPO Op = join(repeat(JWOperator, j), CellOperator);

      // extend the operator to have support over NumCells
      if (NumCells != 0 && j < NumCells-1)
      {
	 Op = join(Op, repeat(identity_mpo(Cell), NumCells-j-1));
      }
      eval.push(element_type(Op));
   }

   UnitCell const& Cell;
   int NumCells;
   std::stack<std::string>& IdentifierStack;
   std::stack<element_type >& eval;
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
   //   static binary_funcs<element_type> binary_funcs_p;

   UnitCellParser(ElementStackType& eval_, 
		  IdentifierStackType& identifier_stack_,
		  UnaryFuncStackType& func_stack_,
		  BinaryFuncStackType& bin_func_stack_,
		  UnitCell const& Cell_, int NumCells_)
      : eval(eval_), identifier_stack(identifier_stack_), 
	func_stack(func_stack_), bin_func_stack(bin_func_stack_), Cell(Cell_), NumCells(NumCells_),
	binary_funcs_p(Cell_) {}
   
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

	 prod_expression = ("prod(" >> expression >> ',' >> expression >> ',' >> quantumnumber >> ')')
	    [push_prod<element_type>(self.Cell, self.identifier_stack, self.eval)];
 
	 bracket_expr = '(' >> expression >> ')';

	 sq_bracket_expr = '[' >> expression >> ']';

	 cell_operator_cell = 
	    identifier 
	    >> (bracket_expr >> (sq_bracket_expr
   [push_local_operator_cell_index(self.Cell, self.NumCells, self.identifier_stack, self.eval)] 
				  |
				  eps_p
   [push_operator_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
				  )
	     | (sq_bracket_expr >> (bracket_expr
   [push_local_operator_index_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
				    |
				    eps_p
   [push_local_operator_no_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)]
				    )
		)
		 )
	    | eps_p[push_operator_no_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)];
	 
	 // cell operator with no cell index specified
	 cell_operator = identifier
	    [push_operator_no_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)];

	 // a cell operator with an index in brackets
	 cell_operator_cell = (identifier >> bracket_expr)
	    [push_operator_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)];

	 // a local operator with no cell index specified
	 local_operator = (identifier >> sq_bracket_expr)
	    [push_local_operator_no_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)];

	 // a local operator with cell index, Op[n](j)
	 local_operator_index_cell = (identifier >> sq_bracket_expr >> bracket_expr)
	    [push_local_operator_index_cell(self.Cell, self.NumCells, self.identifier_stack, self.eval)];

	 // a local operator with cell index, Op(j)[n]
	 local_operator_cell_index = (identifier >> bracket_expr >> sq_bracket_expr)
	    [push_local_operator_cell_index(self.Cell, self.NumCells, self.identifier_stack, self.eval)];
	 
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
							    binary_commutator<element_type> >(self.eval, binary_commutator<element_type>(self.Cell))];
	 
	 factor =
	    imag
	    |   real
	    |   unary_function
	    |   binary_function
	    |   keyword_d[constants_p[push_real<element_type>(self.eval)]]
	    |   prod_expression
	    |   commutator_bracket
	    |   '(' >> expression >> ')'
	    |   ('-' >> factor)[do_negate<element_type>(self.eval)]
	    |   ('+' >> factor)
	    //	    |   local_operator_cell_index
	    //	    |   local_operator_index_cell
	    //	    |   local_operator
	    |   cell_operator_cell
	    //	    |   cell_operator
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
				      binary_multiplication<element_type> >(self.eval, binary_multiplication<element_type>(self.Cell))]
                    |   ('/' >> pow_term)[invoke_binary<element_type, 
					  binary_division<element_type> >(self.eval)]
                    )
	    ;
	 
	 expression =
	    term
	    >> *(  ('+' >> term)[invoke_binary<element_type, 
				 binary_addition<element_type> >(self.eval, binary_addition<element_type>(self.Cell))]
		   |   ('-' >> term)[invoke_binary<element_type, 
				     binary_subtraction<element_type> >(self.eval, binary_subtraction<element_type>(self.Cell))]
		   )
	    ;
      }
      
      rule<ScannerT> expression, term, factor, real, imag, operator_literal, unary_function,
	 binary_function, bracket_expr, quantumnumber, prod_expression, sq_bracket_expr, 
	 cell_operator, cell_operator_cell,
	 local_operator, local_operator_cell_index, local_operator_index_cell, 
	 identifier, pow_term, commutator_bracket;
      rule<ScannerT> const&
      start() const { return expression; }
   };
   
   std::stack<element_type>& eval;
   IdentifierStackType& identifier_stack;
   std::stack<unary_func_type>& func_stack;
   std::stack<binary_func_type>& bin_func_stack;
   UnitCell const& Cell;
   int NumCells;

   binary_funcs<element_type> binary_funcs_p;
};


// global variables (static members of UnitCellParser)
constants UnitCellParser::constants_p;
unary_funcs<UnitCellParser::element_type> UnitCellParser::unary_funcs_p;
//binary_funcs<UnitCellParser::element_type> UnitCellParser::binary_funcs_p;

FiniteMPO
ParseUnitCellOperator(UnitCell const& Cell, int NumCells, std::string const& Str)
{
   typedef UnitCellParser::element_type element_type;

   UnitCellParser::ElementStackType ElementStack;
   UnitCellParser::IdentifierStackType IdentifierStack;
   UnitCellParser::UnaryFuncStackType UnaryFuncStack;
   UnitCellParser::BinaryFuncStackType BinaryFuncStack;

   UnitCellParser Parser(ElementStack, IdentifierStack, UnaryFuncStack, BinaryFuncStack, Cell, NumCells);

   parse_info<> info = parse(Str.c_str(), Parser, space_p);
   if (!info.full)
   {
      PANIC("Operator parser failed, stopped at")(info.stop);
   }

   CHECK(UnaryFuncStack.empty());
   CHECK(BinaryFuncStack.empty());
   CHECK(IdentifierStack.empty());
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
