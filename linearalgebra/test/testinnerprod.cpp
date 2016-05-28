// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testinnerprod.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/vector.h"
#include <complex>
using namespace LinearAlgebra;

typedef std::complex<double> complex;

int main()
{
   Vector<complex> v1 = range(0, 10);
   Vector<complex> v2 = range(20,30);

   TRACE(inner_prod(v1, v2));
   TRACE(inner_prod(herm(v1), v2));
}
