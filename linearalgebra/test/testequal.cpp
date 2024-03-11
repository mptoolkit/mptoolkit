// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testequal.cpp
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
   Matrix<double> A, B;
   Matrix<double> M2(2,2,1);
   Matrix<double> M3(3,3,1);
   Matrix<double> N(2,2,1E-10);
   Matrix<double> Z(2,2,0);

   CHECK(equal(A, B));
   CHECK(equal(M2, M2));
   CHECK(equal(M2+M2, M2*2));
   CHECK(equal(M2,M2+Z));
   CHECK(!equal(M2,M3));
   CHECK(!equal(M2,M2+N));
}
