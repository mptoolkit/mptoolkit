/* -*- C++ -*- $Id: operator-parser.h 1326 2014-03-04 10:56:37Z ianmcc $

  siteoperator-parser.h

  Spirit-based parser for parsing operators.
*/

#if !defined(MPTOOLKIT_LATTICE_SITEOPERATOR_PARSER_H)
#define MPTOOLKIT_LATTICE_SITEOPERATOR_PARSER_H

#include "latticesite.h"

SiteOperator
ParseSiteOperator(LatticeSite const& Site, std::string const& str);

#endif
