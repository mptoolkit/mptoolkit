// -*- C++ -*- $Id$

// mpo operations
// Operations on finite and triangular MPO's

//
// examples of possible syntax:
// sum over a 1-site unit cell (or pretend its a 1-site unit cell)
// sum(j=0, C(j)*C(j+1))
//
// sum over a 2-site unit cell
// sum(j=0 by 2, C(J)*C(j+1) + D(j+1)*D(j+2))
//
// N(k)
// sum(x=0, exp(i*k*x)*CH(x)) * sum(x=0, exp(-i*k*x)*CH(x))
//
// Exponentially decaying interaction
// sum(x=0, sum(y>=0, exp(-y/xi)*CH(x)*C(x+y)))
//
// Finite summation
// sum(x=10 to 20, -(x-15)^2 * N(c))

#if !defined(MPTOOLKIT_MPO_MPO_OPERATIONS_H)
#define MPTOOLKIT_MPO_MPO_OPERATIONS_H

#include "triangular_mpo.h"
#include "finite_mpo.h"

// returns the triangular MPO that corresponds to summing the operator
// over every unit cell.  
// PRECONDITION: Operator.size() is a multiple of UnitCellSize
TriangularMPO
sum_over_unit_cell(FiniteMPO const& Operator, int UnitCellSize);

// Sum the operator over every unit cell, with phase factor
// exp(i*MomentumPerUnitCell) per unit cell.
TriangularMPO
sum_over_unit_cell(FiniteMPO const& Operator, int UnitCellSize,
		   double MomentumPerUnitCell);

// Convert a triangular MPO into a FiniteMPO by restricting the summation
// to a finite size.
// PRECONDITION: Size is a multiple of Operator.size()
// This is sort-of an inverse operation of sum_over_unit_cell.
FiniteMPO
restrict(TriangularMPO const& Operator, int Size);

// Constructs the triangular MPO that represents the infinite sum
// Op1 \otimes Op2
// + Op1 \otimes String \otimes Op2
// + Op1 \otimes String \otimes String \otimes Op2
// + ...
//
// if String.size() < UnitCellSize then the string operator is repeated to make
// a full unit cell.
//
// PRECONDITION: Op1.size() is a multiple of UnitCellSize
// PRECONDITION: Op2.size() is a multiple of UnitCellSize
// PRECONDITION: UnitCellSize is a multiple of String.size()
TriangularMPO
string_product(FiniteMPO const& Op1, FiniteMPO const& String, FiniteMPO const& Op2,
	       int UnitCellSize);

// version of string_product where the string operator is Factor*identity_{UnitCellSize}
TriangularMPO
string_product(FiniteMPO const& Op1, std::complex<double> Factor, FiniteMPO const& Op2,
	       int UnitCellSize);

#endif
