// -*- C++ -*- $Id$

#if !defined(MPTOOLKIT_LATTICE_INFINITELATTICE_PARSER_H)
#define MPTOOLKIT_LATTICE_INFINITELATTICE_PARSER_H

#include "infinitelattice.h"

ProductMPO
ParseProductOperator(InfiniteLattice const& Lattice, std::string const& Str);

std::pair<ProductMPO, InfiniteLattice>
ParseProductOperatorAndLattice(std::string const& Str);

#endif
