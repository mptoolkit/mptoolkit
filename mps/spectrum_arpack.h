// -*- C++ -*- $Id$

// Get the spectrum of an operator using ARPACK

#if !defined(SPECTRUM_H_JZDHCJKYRUYUIFSDHIUO85789289)
#define SPECTRUM_H_JZDHCJKYRUYUIFSDHIUO85789289

#include "infinitewavefunction.h"
#include "packunpack.h"
#include "common/arpackf.h"
#include "tensor/tensor_eigen.h"

LinearAlgebra::Vector<std::complex<double> >
get_spectrum(LinearWavefunction const& Psi, QuantumNumber const& QShift, int NumEigen,
             QuantumNumbers::QuantumNumber const& q, double tol = 1e-10,
             LinearAlgebra::Vector<MatrixOperator>* OutputVectors = NULL, int ncv = 0, int Verbose = 0);

#endif
