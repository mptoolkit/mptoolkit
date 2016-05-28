// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testsliceblas.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER


#define BLAS1_TRACE_DETAILED

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include "common/trace.h"

using namespace LinearAlgebra;

using tracer::typeid_name;

int main()
{
   Vector<double> V(100);
   Vector<double> X(100);

   X[Slice(1,1,1)] = V[Slice(1,1,1)];

   Matrix<double> M(100,100);
   matrix_row(M, 10)[Slice(1,1,1)] = V[Slice(1,1,1)];
}
