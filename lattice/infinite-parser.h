// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/infinite-parser.h
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

#if !defined(MPTOOLKIT_LATTICE_INFINITE_PARSER_H)
#define MPTOOLKIT_LATTICE_INFINITE_PARSER_H

#include "lattice/infinitelattice.h"
#include "mpo/infinite_mpo.h"

InfiniteMPOElement
ParseInfiniteOperator(InfiniteLattice const& Lattice, std::string const& Str,
                      Function::ArgumentList const& Args = Function::ArgumentList());

std::pair<InfiniteMPOElement, InfiniteLattice>
ParseInfiniteOperatorAndLattice(std::string const& Str);

// ProductMPO versions

ProductMPO
ParseProductOperator(InfiniteLattice const& Lattice, std::string const& Str,
                      Function::ArgumentList const& Args = Function::ArgumentList());


std::pair<ProductMPO, InfiniteLattice>
ParseProductOperatorAndLattice(std::string const& Str);

// TriangularMPO versions

TriangularMPO
ParseTriangularOperator(InfiniteLattice const& Lattice, std::string const& Str,
                        Function::ArgumentList const& Args = Function::ArgumentList());


std::pair<TriangularMPO, InfiniteLattice>
ParseTriangularOperatorAndLattice(std::string const& Str);

// complex versions

std::complex<double>
ParseInfiniteNumber(InfiniteLattice const& Lattice, std::string const& Str,
                    Function::ArgumentList const& Args = Function::ArgumentList());

namespace Parser
{

inline
std::complex<double>
as_number(InfiniteMPO const& x)
{
   return x.as_complex();
}

} // namespace Parser

using namespace Parser;

#endif
