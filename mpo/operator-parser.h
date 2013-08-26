/* -*- C++ -*- $Id$

  operator-parser.h

  Spirit-based parser for parsing operators.
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
