// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/unitcell-parser.h
//
// Copyright (C) 2014-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "common/sha256.h"

extern InfiniteLattice const* ILattice;

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
