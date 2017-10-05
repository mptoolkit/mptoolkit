// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testeigensymmetric.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
#include "linearalgebra/eigen.h"

#include <iostream>
#include <iomanip>

using namespace std;
using namespace LinearAlgebra;

int main()
{
  Matrix<double> A(2,2, 0);
  Matrix<double> B(2,2, 0);
  B(0,0) = 4;  B(1,1) = 4;  B(1,0) = B(0,1) = 1;

  A(0,0) = 1; A(1,1) = 1;
  A(1,0) = A(0,1) = 4;

  Vector<double> Evalues;
  Matrix<double> Evectors;

  GeneralizedEigenHermitian(A, B, Evalues, Evectors, Range(0,2));

  Matrix<double> ExpectedEvectors(2,2);
  ExpectedEvectors(0,0) = -1.0 / std::sqrt(6.0);
  ExpectedEvectors(0,1) = 1.0 / std::sqrt(6.0);
  ExpectedEvectors(1,0) = ExpectedEvectors(1,1) = 1.0 / std::sqrt(10.0);
  Vector<double> ExpectedEvalues(2);
  ExpectedEvalues[0] = -1;
  ExpectedEvalues[1] = 1;

  CHECK(equal(Evalues, ExpectedEvalues));
  CHECK(equal(Evectors, ExpectedEvectors));

}
