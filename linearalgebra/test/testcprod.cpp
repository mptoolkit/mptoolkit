// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testcprod.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
  Matrix<std::complex<double> > M(3,4,std::complex<double>(2.0,1.0));
  Matrix<std::complex<double> > N(3,4,4.0);
  Matrix<std::complex<double> > R; R = M*transpose(conj(N));
  std::cout << R;
}
