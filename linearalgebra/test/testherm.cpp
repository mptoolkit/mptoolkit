// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testherm.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/matrix.h"
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/matrix_utility.h"
#include <complex>

using namespace LinearAlgebra;

int main()
{
   typedef Matrix<std::complex<double> > T;
   T M1 = random_matrix<T::value_type>(2,2);
   SparseMatrix<T> M(4,5);
   M(0,2) = M1;
   TRACE(conj(-M));
   TRACE(herm(M));
   TRACE(transform(M, Herm<T>()));
   SparseMatrix<T> N(herm(M));
}
