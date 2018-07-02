// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix_utility.h
//
// Copyright (C) 2004-2018 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_BLAS_MATRIX_UTILITY_H)
#define MPTOOLKIT_BLAS_MATRIX_UTILITY_H

#include "matrix.h"

namespace blas
{

// generates a random matrix with elements uniformly distributed in the range (-1,1)
template <typename T>
Matrix<T> random_matrix(int r, int c);

} // namespace blas

#include "matrix_utility.icc"

#endif
