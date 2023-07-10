// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/momentum_operations.h
//
// Copyright (C) 2015-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "mpo/basic_triangular_mpo.h"
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

// scale the matrices by a factor
void
ScalePoly(KMatrixPolyType& K, std::complex<double> x);

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
             BasicTriangularMPO const& Op,
             LinearWavefunction const& Psi,
             QuantumNumber const& QShift, int Column);

// Calculate the polynomial of overlaps of the E matrix with some operator (typically the density matrix)
//  |---|
//  E* Rho
//  |---|
ComplexPolyType
ExtractOverlap(MatrixPolyType const& E, MatrixOperator const& Rho);

// Extract overlaps, decomposed into momenta
KComplexPolyType
ExtractOverlap(KMatrixPolyType const& E, MatrixOperator const& Rho);

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
add_triple_prod(MatrixPolyType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                MatrixPolyType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp);

void
add_triple_prod(KMatrixPolyType& Result, std::complex<double> Factor,
                HermitianProxy<MatrixOperator> const& x,
                KMatrixPolyType const& E,
                MatrixOperator const& y,
                QuantumNumber const& qxy,
                QuantumNumber const& qEp);

void
add_triple_prod(MatrixPolyType& Result, std::complex<double> Factor,
            MatrixOperator const& x,
            MatrixPolyType const& E,
            HermitianProxy<MatrixOperator> const& y,
            QuantumNumber const& qxy,
            QuantumNumber const& qEp);

void
add_triple_prod(KMatrixPolyType& Result, std::complex<double> Factor,
            MatrixOperator const& x,
            KMatrixPolyType const& E,
            HermitianProxy<MatrixOperator> const& y,
            QuantumNumber const& qxy,
            QuantumNumber const& qEp);

std::vector<KMatrixPolyType>
contract_from_left(OperatorComponent const& M,
                   HermitianProxy<StateComponent> const& A,
                   std::vector<KMatrixPolyType> const& E,
                   StateComponent const& B);

template <typename EType>
std::vector<EType>
contract_from_left_mask(OperatorComponent const& M,
                        HermitianProxy<StateComponent> const& A,
                        std::vector<EType> const& E,
                        StateComponent const& B,
                        std::vector<int> const& OutMask,
                        std::vector<int> const& InMask);

template <typename EType>
std::vector<EType>
contract_from_right_mask(HermitianProxy<OperatorComponent> const& M,
                        StateComponent const& A,
                        std::vector<EType> const& F,
                        HermitianProxy<StateComponent> const& B,
                        std::vector<int> const& OutMask,
                        std::vector<int> const& InMask);

std::vector<KMatrixPolyType>
inject_left(std::vector<KMatrixPolyType> const& In,
            LinearWavefunction const& Psi1,
            GenericMPO const& Op,
            LinearWavefunction const& Psi2);

std::vector<MatrixPolyType>
inject_left_mask(std::vector<MatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int>> const& Mask);

std::vector<MatrixPolyType>
inject_right_mask(std::vector<MatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int>> const& Mask);

std::vector<KMatrixPolyType>
inject_left_mask(std::vector<KMatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int>> const& Mask);

std::vector<KMatrixPolyType>
inject_right_mask(std::vector<KMatrixPolyType> const& In,
                 LinearWavefunction const& Psi1,
                 QuantumNumber const& QShift,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int>> const& Mask);

//
// implementation
//

template <typename EType>
std::vector<EType>
contract_from_left_mask(OperatorComponent const& M,
                        HermitianProxy<StateComponent> const& A,
                        std::vector<EType> const& E,
                        StateComponent const& B,
                        std::vector<int> const& OutMask,
                        std::vector<int> const& InMask)
{
  std::vector<EType> Result(M.Basis2().size());

  // Iterate over the components in M, first index
  for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M); I; ++I)
  {
     // skip over masked components
     if (!InMask[I.index()])
        continue;

     // second index in M
     for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
     {
        // skip over masked components
        if (!OutMask[J.index2()])
           continue;

        // Iterate over the irreducible components of M(I,J)
        for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
        {
           // *k is an irreducible operator.  Iterate over the components of this operator
           for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
           {
              for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type
                      S = iterate(R); S; ++S)
              {
                 add_triple_prod(Result[J.index2()], *S,
                                 herm(A.base()[S.index1()]),
                                 E[J.index1()],
                                 B[S.index2()],
                                 k->TransformsAs(),
                                 M.Basis2()[J.index2()]);
              }
           }
        }
     }
  }
  return Result;
}

template <typename FType>
std::vector<FType>
contract_from_right_mask(HermitianProxy<OperatorComponent> const& M,
                        StateComponent const& A,
                        std::vector<FType> const& F,
                        HermitianProxy<StateComponent> const& B,
                        std::vector<int> const& OutMask,
                        std::vector<int> const& InMask)
{
   PRECONDITION_EQUAL(M.base().LocalBasis1(), A.LocalBasis());
   PRECONDITION_EQUAL(M.base().LocalBasis2(), B.base().LocalBasis());

  std::vector<FType> Result(M.base().Basis1().size());

   // Iterate over the components in M, first index
   for (LinearAlgebra::const_iterator<OperatorComponent>::type I = iterate(M.base()); I; ++I)
   {
      // skip over masked components
      if (!OutMask[I.index()])
         continue;

      // second index in M
      for (LinearAlgebra::const_inner_iterator<OperatorComponent>::type J = iterate(I); J; ++J)
      {
         // skip over masked components
         if (!InMask[J.index2()])
            continue;

         // Iterate over the irreducible components of M(I,J)
         for (SimpleRedOperator::const_iterator k = J->begin(); k != J->end(); ++k)
         {
            // *k is an irreducible operator.  Iterate over the components of this operator
            for (LinearAlgebra::const_iterator<SimpleOperator>::type R = iterate(*k); R; ++R)
            {
               for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type
                       S = iterate(R); S; ++S)
               {
                  add_triple_prod(Result[J.index1()], herm(*S),
                                  A[S.index1()],
                                  F[J.index2()],
                                  herm(B.base()[S.index2()]),
                                  k->TransformsAs(),
                                  M.base().Basis1()[J.index1()]);
               }
            }
         }
      }
   }
   return Result;
}

#endif
