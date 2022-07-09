// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/write-matrix.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// Functions do write matrices, operators and MPS in various formats.
// Currently supported: python, MATLAB

#if !defined(MPTOOLKIT_MP_ALGORITHMS_WRITE_MATRIX_H)
#define MPTOOLKIT_MP_ALGORITHMS_WRITE_MATRIX_H

#include "wavefunction/mpwavefunction.h"
#include <iostream>

void WriteMatrixFormatMATLAB(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M, bool Quiet = false);

void WriteRealMatrixFormatMATLAB(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M, bool quiet = false);

void WriteMPS_MATLAB(std::ostream& out, LinearWavefunction const& Psi, MatrixOperator const& Rho, bool Quiet = false);

void WriteMatrixFormatPython(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                            std::string Prefix = "", bool Quiet = false);

void WriteRealMatrixFormatPython(std::ostream& out, LinearAlgebra::Matrix<std::complex<double>> const& M,
                                 std::string Prefix = "", bool quiet = false);

void WriteMPS_Python(std::ostream& out, LinearWavefunction const& Psi, MatrixOperator const& Rho, bool Quiet = false);

#endif // defined(MPTOOLKIT_MP_ALGORITHMS_WRITE_MATRIX_H)
