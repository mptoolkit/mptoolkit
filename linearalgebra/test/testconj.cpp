// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testconj.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using namespace std;

int main()
{
  double y = 3;
  double x = conj(y);
  CHECK_EQUAL(x,y);

  Vector<double> MV(3,-1);
  Vector<double> NV = conj(MV);
  CHECK_CLOSE(MV, NV);

  Vector<complex<double> > ZV(3,complex<double>(0,1));
  Vector<complex<double> > CV(conj(ZV));
  CHECK_CLOSE(CV, -ZV);

  Matrix<double> M(3,3,-1);
  Matrix<double> N = conj(M);
  CHECK_CLOSE(M, N);

  Matrix<complex<double> > Z(3,3,complex<double>(0,1));
  Matrix<complex<double> > C(conj(Z));
  CHECK_CLOSE(C, -Z);
}
