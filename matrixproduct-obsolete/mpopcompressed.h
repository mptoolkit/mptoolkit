// -*- C++ -*- $Id$
//
// Run-length compressed matrix product operator,
// and associated low-level functions
//

#if !defined(MPOPCOMPRESSED_H_JDCHJKEHY589758YUER89H489)
#define MPOPCOMPRESSED_H_JDCHJKEHY589758YUER89H489

#include "common/runlengthcompressed.h"
#include "pstream/pstream.h"
#include "quantumnumbers/quantumnumber.h"
#include "mpopcomponent.h"

typedef run_length_compressed<MPOpComponent> MPOpCompressed;
typedef run_length_repeat<MPOpComponent> MPOpRepeat;
typedef run_length_array<MPOpComponent> MPOpArray;

// Splits a compressed operator into two pieces, the first piece has size
// Loc, the second piece has size x.size()-Loc.
inline
std::pair<MPOpCompressed, MPOpCompressed>
split_operator(MPOpCompressed const& x, int Loc)
{
   return split(x, Loc);
}

// Calculates the product of operators, from left to right.  Here,
// B is a ProductBasis for a.front().Basis1() x b.front().Basis2(),
// and C is some scalar operator that acts from left to right.
// On exit, C is shifted to the far right hand side, and B is the ProductBasis
// for the original right-most Basis2() of a and b.
MPOpCompressed do_prod(SimpleOperator& C, ProductBasis<BasisList, BasisList>& B,
                       MPOpCompressed const& a, MPOpCompressed const& b);

// Calculates the sum of operators, from left to right.
// Initial B is a SumBasis for a.front().Basis1() + b.front().Basis1(),
// and C is an operator with C.Basis2() == B.
// On exit, C is shifted to the far right hand side, and B is the SumBasis
// for the original right-most Basis2() of a and b.
MPOpCompressed do_sum(SimpleOperator& C, SumBasis<BasisList>& B,
                      MPOpCompressed const& a, MPOpCompressed const& b);

// TODO: The implementations of do_prod and do_sum are almost identical.
// Can we simplify somehow?

// Injects the operator C into a, from left to right.  On exit,
// C is shifted to the right hand side.
MPOpCompressed inject_left(SimpleOperator& C, MPOpCompressed const& a);

// Injects the operator C into a, from right to left.  On exit,
// C is shifted to the left hand side.
MPOpCompressed inject_right(MPOpCompressed const& a, SimpleOperator& C);

// Multiplies the left-most component of a by C.
MPOpCompressed multiply_left(SimpleOperator const& C, MPOpCompressed const& a);

// Projects onto a component n1 of the left-most basis and n2 of the right-most basis
MPOpCompressed project_component(MPOpCompressed const& Res, int n1, int n2);

// The adjoint of an operator
MPOpCompressed adjoint(MPOpCompressed const& x);

// The conjugate of an operator
MPOpCompressed conj(MPOpCompressed const& x);

std::ostream& operator<<(std::ostream& out, MPOpCompressed const& x);

// some helper functions, dont really belong here but are used in the
// implementation here, and in mpoperator.cpp

// helper function to construct an operator that merges repeated
// quantum numbers in the basis.  
SimpleOperator CollapseBasis(BasisList const& b);

SimpleOperator RemoveEmptyRows(SimpleOperator const& c);

#endif
