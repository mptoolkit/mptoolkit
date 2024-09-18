// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/spectrum_arpack.h
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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
