/* -*- C++ -*- $Id: operator-parser.h 1326 2014-03-04 10:56:37Z ianmcc $

  siteoperator-parser.h

  Spirit-based parser for parsing operators.
*/

#if !defined(MPTOOLKIT_LATTICE_SITEOPERATOR_PARSER_H)
#define MPTOOLKIT_LATTICE_SITEOPERATOR_PARSER_H

#include "latticesite.h"
#include "function.h"
#include <boost/variant.hpp>

typedef boost::variant<SiteOperator, std::complex<double> >
SiteElementType;

// Parse a string, resulting in either an operator or a c-number
SiteElementType
ParseSiteElement(LatticeSite const& Site, 
		 std::string const& Str,
		 Function::ArgumentList const& Args = Function::ArgumentList());

// Parse a string, where we expect the result to be a
// SiteOperator
SiteOperator
ParseSiteOperator(LatticeSite const& Site, 
		  std::string const& Str,
		  Function::ArgumentList const& Args = Function::ArgumentList());

// Parse a string, where we expect the result to be a
// c-number
std::complex<double>
ParseSiteNumber(LatticeSite const& Site, 
		std::string const& Str,
		Function::ArgumentList const& Args = Function::ArgumentList());

namespace Parser
{

// conversion of an element_type to a c-number
inline
std::complex<double>
as_number(boost::variant<SiteOperator, std::complex<double> > const& x)
{
   return boost::get<std::complex<double> >(x);
}

} // namespace Parser

#endif
