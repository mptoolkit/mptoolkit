// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/spectrum_arpack.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

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
