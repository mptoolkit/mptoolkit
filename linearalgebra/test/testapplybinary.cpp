// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testapplybinary.cpp
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
  Matrix<double> M(3,4,2.0);
  Matrix<double> N(3,4,4.0);
  M(0,1) = 1;

  Matrix<double> MNCheck(3,4,8.0);
  MNCheck(0,1) = 4;

  Matrix<double> MNElementProd = apply_func(BinaryOperator<Multiplication, double, double>(), M, N);

  CHECK_EQUAL(MNElementProd, MNCheck);
}
