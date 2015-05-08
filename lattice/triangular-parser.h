// -*- C++ -*- $Id$

#if !defined(MPTOOLKIT_LATTICE_INFINITELATTICE_PARSER_H)
#define MPTOOLKIT_LATTICE_INFINITELATTICE_PARSER_H

#include "infinitelattice.h"

TriangularMPO
ParseTriangularOperator(InfiniteLattice const& Lattice, std::string const& Str);

std::pair<TriangularMPO, InfiniteLattice>
ParseTriangularOperatorAndLattice(std::string const& Str);

#endif
