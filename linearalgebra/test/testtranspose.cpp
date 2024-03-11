// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testtranspose.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "linearalgebra/matrix.h"
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
  CHECK_EQUAL(MV,NV);

  Vector<complex<double> > ZV(3,complex<double>(0,1));
  Vector<complex<double> > CV(conj(ZV));
  CHECK_EQUAL(CV, -ZV);

  Matrix<double> M(3,3,-1);
  Matrix<double> N = conj(M);
  CHECK_EQUAL(M,N);

  Matrix<complex<double> > Z(3,3,complex<double>(0,1));
  Matrix<complex<double> > C(conj(Z));
  CHECK_EQUAL(C, -Z);
}
