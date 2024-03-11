// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/takagi.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MPTOOLKIT_LINEARALGEBRA_TAKAGI_H)
#define MPTOOLKIT_LINEARALGEBRA_TAKAGI_H

#include "linearalgebra/matrix.h"
#include "linearalgebra/diagonalmatrix.h"

//   TakagiFactor factorizes a complex symmetric n-by-n matrix
//   Input:	n, A = n-by-n matrix, complex symmetric
//		(only the upper triangle of A needs to be filled),
//   Output:	d = vector of diagonal values,
//		U = transformation matrix, unitary (U^-1 = U^+),
//   these fulfill
//	d = U^* A U^+,  A = U^T d U,  U^* A = d U

std::tuple<LinearAlgebra::Matrix<std::complex<double>>, LinearAlgebra::DiagonalMatrix<double>>
TakagiFactor(LinearAlgebra::Matrix<std::complex<double>> A);

#endif
