// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/matrix-lapack.h
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

#if !defined(MPTOOLKIT_BLAS_MATRIX_LAPACK_H)
#define MPTOOLKIT_BLAS_MATRIX_LAPACK_H

#include <complex>

//
// LAPACK wrappers for ordinary matrices
//

namespace blas
{

void DiagonalizeSymmetric(int Size, double* Data, int LeadingDim, double* Eigen);

void DiagonalizeHermitian(int Size, std::complex<double>* Data, int LeadingDim, double* Eigen);

void SingularValueDecomposition(int Rows, int Cols, double* Data, int LeadingDim, double* Dvec,
				double* Umat, int ldU, double* Vmat, int ldV);

} // namespace blas


#endif
