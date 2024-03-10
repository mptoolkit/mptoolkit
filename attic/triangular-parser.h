// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/triangular-parser.h
//
// Copyright (C) 2015-2017 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
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
