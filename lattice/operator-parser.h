// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/operator-parser.h
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
/* -*- C++ -*- $Id$

  operator-parser.h

  Spirit-based parser for parsing operators.

  *** OBSOLETE ***
*/

#if !defined(OPERATOR_PARSER_H_JDSHCIUYY876Y76Y743YT66IUREYT473YTRKH)
#define OPERATOR_PARSER_H_JDSHCIUYY876Y76Y743YT66IUREYT473YTRKH

#include "matrixproduct/mpoperatorlist.h"
#include <boost/algorithm/string.hpp>
#include <utility>

//
// ParseOperator
//
// Parse an expression of the form lattice:operator
// and return the corresponding MPOperator.

std::pair<OperatorList, MPOperator>
ParseLatticeAndOperator(std::string const& str);

inline
MPOperator ParseOperator(std::string const& str)
{
   return ParseLatticeAndOperator(str).second;
}

#endif
