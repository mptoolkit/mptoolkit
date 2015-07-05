// -*- C++ -*- $Id: product-parser.h 1461 2015-05-08 17:01:39Z ianmcc $

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
