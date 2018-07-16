// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// blas/exponential.h
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

//
// A simple dense matrix class designed for scalar types.
//

#if !defined(MPTOOLKIT_BLAS_EXPONENTIAL_H)
#define MPTOOLKIT_BLAS_EXPONENTIAL_H

#include "matrix.h"

namespace blas
{

template <typename T, typename Tag>
Matrix<T, Tag>
exp(Matrix<T, Tag> const& m);

Matrix<float>
exp(Matrix<float> const& m);

Matrix<std::complex<float>>
exp(Matrix<std::complex<float>> const& m);

Matrix<double>
exp(Matrix<double> const& m);

Matrix<std::complex<double>>
exp(Matrix<std::complex<double>> const& m);

} // namepsace blas

#include "exponential.icc"

#endif
