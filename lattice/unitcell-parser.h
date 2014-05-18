/* -*- C++ -*- $Id$

  unitcell-parser.h

  Spirit-based parser for parsing operators.
*/

#if !defined(MPTOOLKIT_LATTICE_UNITCELL_PARSER_H)
#define MPTOOLKIT_LATTICE_UNIT_PARSER_H

#include "unitcell.h"

FiniteMPO
ParseUnitCellOperator(UnitCell const& Cell, std::string const& str);

#endif
