// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testtrace.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
  Matrix<double> A(2,2);
  A(0,0) = 3;
  A(1,1) = -1;
  A(1,0) = 2;
  A(0,1) = 8;
  CHECK_EQUAL(trace(A), 2);
  CHECK_EQUAL(inner_prod(A,A), 78);
  CHECK_EQUAL(inner_prod(A,trans(A)), 42);
}
