// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mpo/operator_utilities.h
//
// Copyright (C) 2016 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_MPO_OPERATOR_UTILITIES_H)
#define MPTOOLKIT_MPO_OPERATOR_UTILITIES_H

#include "tensor/tensor.h"
#include "tensor/reducible.h"
#include "linearalgebra/vector.h"
#include "mps/state_component.h"

// default epsilon for detecting whether an eigenvalue is equal to 1 for
// operator classifications
double const DefaultClassifyUnityEpsilon = 1E-14;

// a simple reducible operator typedef
typedef ReducibleTensor<std::complex<double>, BasisList, BasisList> SimpleRedOperator;

int
linear_dimension(SimpleOperator const& c);

LinearAlgebra::Vector<std::complex<double>>
linearize(SimpleOperator const& c);

// Linearize a SimpleRed operator into a vector, of full dimensions (ie, every degree of freedom is
// represented), in the given quantum number sectors.
LinearAlgebra::Vector<std::complex<double>>
linearize(SimpleRedOperator const& c, QuantumNumberList const& QList);

int
linear_dimension(SimpleRedOperator const& c, QuantumNumberList const& QList);

typedef IrredTensor<std::complex<double> > SimpleOperator;

typedef IrredTensor<double> RealSimpleOperator;

typedef IrredTensor<LinearAlgebra::Matrix<std::complex<double> >,
                            VectorBasis,
                            VectorBasis> MatrixOperator;

// Constructs the swap gate, that maps from the tensor product basis of B1 * B2 into
// the basis B2 * B1, such that state |i,j> maps into |j,i>
SimpleOperator
swap_gate(BasisList const& B1, BasisList const& B2,
          ProductBasis<BasisList, BasisList> const& Basis_21,
          ProductBasis<BasisList, BasisList> const& Basis_12);



inline
SimpleOperator
swap_gate(BasisList const& B1, BasisList const& B2)
{
   return swap_gate(B1, B2, ProductBasis<BasisList, BasisList>(B2,B1), ProductBasis<BasisList, BasisList>(B1,B2));
}

// Fermionic swap gate that inserts (-1)^{N_1 * N_2} into the swap gate.
// Parity1 and Parity2 are the diagonal components of the P operator, ie Parity[i] == 1
// iff B1(i) is bosonic, and Parity1[i] == -1 iff B1(i) is fermionic.
SimpleOperator
swap_gate_fermion(BasisList const& B1, LinearAlgebra::Vector<double> const& Parity1,
                  BasisList const& B2, LinearAlgebra::Vector<double> const& Parity2,
                  ProductBasis<BasisList, BasisList> const& Basis_21,
                  ProductBasis<BasisList, BasisList> const& Basis_12);

inline
SimpleOperator
swap_gate(BasisList const& B1, LinearAlgebra::Vector<double> const& Parity1,
          BasisList const& B2, LinearAlgebra::Vector<double> const& Parity2)
{
   return swap_gate(B1, B2, ProductBasis<BasisList, BasisList>(B2,B1), ProductBasis<BasisList, BasisList>(B1,B2));
}

// helper function to construct an operator that merges repeated
// quantum numbers in the basis.
SimpleOperator CollapseBasis(BasisList const& b);

// helper function to construct an operator that projects onto a single quantum number component
SimpleOperator ProjectBasis(BasisList const& b, QuantumNumbers::QuantumNumber const& q);

// Utility function - if X is proportional to identity operator then return the constant of
// proportionality.  Otherwise return 0.0
std::complex<double> PropIdent(SimpleOperator const& X, double UnityEpsilon);

inline
std::complex<double> PropIdent(SimpleOperator const& X)
{
   return PropIdent(X, DefaultClassifyUnityEpsilon);
}

#endif
