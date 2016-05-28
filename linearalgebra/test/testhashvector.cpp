// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testhashvector.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/hashvector.h"
#include "testvectorsparse.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrixmatrixmultiplication.h"

int main()
{
   test_real_sparse_vector<LinearAlgebra::HashVector<double> >();

#if 0
   LinearAlgebra::HashVector<LinearAlgebra::Matrix<double> > v1(10), v2(10);
   TRACE(typeid(parallel_prod(v1,v2)).name());
#endif
}
