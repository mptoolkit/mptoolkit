// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// attic/triangular-parser.h
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

#if !defined(MPTOOLKIT_LATTICE_INFINITELATTICE_PARSER_H)
#define MPTOOLKIT_LATTICE_INFINITELATTICE_PARSER_H

#include "infinitelattice.h"

BasicTriangularMPO
ParseTriangularOperator(InfiniteLattice const& Lattice, std::string const& Str);

std::pair<BasicTriangularMPO, InfiniteLattice>
ParseTriangularOperatorAndLattice(std::string const& Str);

#endif
