/* -*- C++ -*-

  unitcell-parser.h

  Spirit-based parser for parsing operators.

  Parses an operator that has finite support over some set of unit cells.

  If the operator has support over only one unit cell, then the
  cell index can be omitted (it is implicitly zero).

  If the NumCells parameter is zero, then the operator can have support over an arbitrary number
  of unit cells; the MPO will be extended to suit.
*/

#if !defined(MPTOOLKIT_LATTICE_UNITCELL_PARSER_H)
#define MPTOOLKIT_LATTICE_UNITCELL_PARSER_H

#include "lattice/unitcell.h"
#include "lattice/unitcell_mpo.h"
#include "lattice/infinitelattice.h"

typedef boost::variant<UnitCellMPO, std::complex<double> >
UnitCellElementType;

UnitCellElementType
ParseUnitCellElement(UnitCell const& Cell, int NumCells, std::string const& str,
		     Function::ArgumentList const& Args = Function::ArgumentList());

UnitCellMPO
ParseUnitCellOperator(UnitCell const& Cell, int NumCells, std::string const& str,
		      Function::ArgumentList const& Args = Function::ArgumentList());

std::complex<double>
ParseUnitCellNumber(UnitCell const& Cell, int NumCells, std::string const& str,
		    Function::ArgumentList const& Args = Function::ArgumentList());

std::pair<UnitCellMPO, InfiniteLattice>
ParseUnitCellOperatorAndLattice(std::string const& Str);
   
namespace Parser
{

// conversion of an element_type to a c-number
inline
std::complex<double>
as_number(boost::variant<std::complex<double>, UnitCellMPO> const& x)
{
   return boost::get<std::complex<double> >(x);
}

} // namespace Parser



#endif
