// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testswap.cpp
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

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> V(2);
   V[0] = 3;
   V[1] = 4;

   Vector<double> V2(3);
   V2[0] = -3;
   V2[1] = -4;
   V2[2] = 5;

   Vector<double> U = V;
   Vector<double> U2 = V2;

   swap(V, V2);
   CHECK_EQUAL(V, U2);
   CHECK_EQUAL(V2, U);
}
