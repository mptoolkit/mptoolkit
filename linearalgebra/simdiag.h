// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/simdiag.h
//
// Copyright (C) 2012-2022 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_LINEARALGEBRA_SIMDIAG_H)
#define MPTOOLKIT_LINEARALGEBRA_SIMDIAG_H

#include "linearalgebra/matrix.h"
#include <complex>

// Algorithm to simultaneously diagonalize a set of Hermitian matrices

std::tuple<LinearAlgebra::Matrix<std::complex<double>>, std::vector<LinearAlgebra::Vector<double>>>
SimultaneousDiagonalizeHermitian(std::vector<LinearAlgebra::Matrix<std::complex<double>>> M);

#endif
