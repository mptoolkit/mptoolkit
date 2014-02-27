// -*- C++ -*- $Id$
//
// functions for operators at finite momentum

#if !defined(MPTOOLKIT_MPS_MOMENTUM_OPERATIONS_H)
#define MPTOOLKIT_MPS_MOMENTUM_OPERATIONS_H

#include "common/polynomial.h"
#include "mps/linearwavefunction.h"
#include "mpo/generic_mpo.h"
#include "mpo/triangular_mpo.h"

//
// Polynomial operations, for triangular expectation values (no momentum)
//

// polynomial with matrix coefficients
typedef Polynomial<MatrixOperator> MatrixPolyType;

// polynomial with complex coefficients
typedef Polynomial<std::complex<double> > ComplexPolyType;

// Comparitor for complex numbers.  This is so that we can put them in a map,
// the choice of comparison operation is arbitrary
struct CompareComplex
{
   typedef std::complex<double> first_argument_type;
   typedef std::complex<double> second_argument_type;
   typedef bool result_type;
   bool operator()(std::complex<double> const& x, std::complex<double> const& y) const
   {
      return (x.real() < y.real()) || (x.real() == y.real() && x.imag() < y.imag());
   }
};

// Momentum-dependent complex polynomial
typedef std::map<std::complex<double>, ComplexPolyType, CompareComplex> KComplexPolyType;

// momentum-dependent matrix polynomial,
// this represents an E matrix
typedef std::map<std::complex<double>, MatrixPolyType, CompareComplex> KMatrixPolyType;


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
operator_prod(HermitianProxy<OperatorComponent> const& M,
              HermitianProxy<StateComponent> const& A,
              std::vector<KMatrixPolyType> const& E, 
              StateComponent const& B);

std::vector<KMatrixPolyType>
operator_prod(HermitianProxy<OperatorComponent> const& M,
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
