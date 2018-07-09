// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-blas.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

// standard BLAS functions

#if !defined(MPTOOLKIT_BLAS_MATRIX_BLAS_H)
#define MPTOOLKIT_BLAS_MATRIX_BLAS_H

// generic template versions

#include "matrix-blas-generic.h"

// optimized versions for float/double/complex

// determine which BLAS library we use
#if defined(HAVE_OPENBLAS)
// openblas
#include "matrix-openblas.h"
#elif defined(HAVE_MKL)
// MKL
#include "matrix-mkl.h"

// last possibility is reference BLAS
#eif !defined(HAVE_BLAS)
#error "cannot identify a BLAS library!"
#else
// reference BLAS
#include "matrix-blasreference.h"
#endif

#endif
