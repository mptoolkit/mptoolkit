// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testapplybinary.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/matrix.h"
#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
using namespace std;

int main()
{
  Matrix<double> M(3,4,2.0);
  Matrix<double> N(3,4,4.0);
  M(0,1) = 1;

  Matrix<double> MNCheck(3,4,8.0);
  MNCheck(0,1) = 4;

  Matrix<double> MNElementProd = apply_func(BinaryOperator<Multiplication, double, double>(), M, N);

  CHECK_EQUAL(MNElementProd, MNCheck);
}
