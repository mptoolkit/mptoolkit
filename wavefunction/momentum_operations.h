// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/momentum_operations.h
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
//
// functions for operators at finite momentum

#if !defined(MPTOOLKIT_MPS_MOMENTUM_OPERATIONS_H)
#define MPTOOLKIT_MPS_MOMENTUM_OPERATIONS_H

#include "common/polynomial.h"
#include "wavefunction/linearwavefunction.h"
#include "mpo/generic_mpo.h"
#include "mpo/triangular_mpo.h"
#include "common/angle_map.h"

//
// Polynomial operations, for triangular expectation values (no momentum)
//

// polynomial with matrix coefficients
typedef Polynomial<MatrixOperator> MatrixPolyType;

// polynomial with complex coefficients
typedef Polynomial<std::complex<double> > ComplexPolyType;


// Momentum-dependent complex polynomial.
typedef angle_map<ComplexPolyType> KComplexPolyType;

// momentum-dependent matrix polynomial,
// this represents an E matrix
typedef angle_map<MatrixPolyType> KMatrixPolyType;


MatrixPolyType
delta_shift(MatrixPolyType const& In, QuantumNumber const& QShift);

MatrixPolyType
inject_left(MatrixPolyType const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2);

std::vector<MatrixPolyType>
inject_left(std::vector<MatrixPolyType> const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2);

// Calculates result' = C[column] = sum_{j > Column} E[j] * T_Op(j, Column)
// Assumes that E[j] is defined, for j > Column
MatrixPolyType
MultiplyLeft(std::vector<MatrixPolyType> const& E,
             TriangularMPO const& Op,
             LinearWavefunction const& Psi,
             QuantumNumber const& QShift, int Column);

// Calculate the polynomial of overlaps of the E matrix with some operator (typically the density matrix)
//  |---|
//  E* Rho
//  |---|
ComplexPolyType
ExtractOverlap(MatrixPolyType const& E, MatrixOperator const& Rho);

//
// With momentum
//
// For finite momentum, we extend the MatrixPolyType to be a map
// from std::complex<double> to polynomials.
// An alternative way would be to store the phase angle [0,2*pi), although
// that seems to give no advantage?
//

// delta-shift all components of a MomentumPolynomial operator
KMatrixPolyType
delta_shift(KMatrixPolyType const& In, QuantumNumber const& QShift);

std::vector<KMatrixPolyType>
delta_shift(std::vector<KMatrixPolyType> const& In, QuantumNumber const& QShift);

void
add_triple_prod(KMatrixPolyType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                KMatrixPolyType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp);

std::vector<KMatrixPolyType>
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   std::vector<KMatrixPolyType> const& E,
                   StateComponent const& B);

std::vector<KMatrixPolyType>
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   std::vector<KMatrixPolyType> const& E,
                   StateComponent const& B,
                   std::vector<int> const& OutMask,
                   std::vector<int> const& InMask);

std::vector<KMatrixPolyType>
inject_left(std::vector<KMatrixPolyType> const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2);

std::vector<KMatrixPolyType>
inject_left_mask(std::vector<KMatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask);

#endif
