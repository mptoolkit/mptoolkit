// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testnorminf.cpp
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

#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
  Vector<double> v(3);
  v[0] = -10;
  v[1] = 3;
  v[2] = 7;
  CHECK_EQUAL(norm_inf(v), 10);

  Vector<std::complex<double> > vc(10, std::complex<double>(3,4));
  CHECK_EQUAL(norm_inf(vc), 5);
}
