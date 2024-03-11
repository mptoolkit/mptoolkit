// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testindex.cpp
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
#include "linearalgebra/fixedvector.h"

#include <iostream>

using namespace LinearAlgebra;

int main()
{

   Vector<double> V(10);    // vector of length 10, uninitialized values
   Vector<double> U(10, 1.0); // vector length 10, initialized to 1.0

   fill(U,1.0);  // an alternative way to set all elements to 1.0

   TRACE(U);       // here, U = (1,1,1,1,1,1,1,1,1,1)

   Vector<int> FirstHalf  = Range(0, 5);        // range [0, 5), ie. elements 0,1,2,3,4
   Vector<int> SecondHalf = Range(5, 10);       // range [5, 10), ie. elements 5,6,7,8,9

   Vector<int> EvenElements = Slice(0,5,2);     // slice starting at element 0, length 5, stride 2
   Vector<int> OddElements = Slice(1,5,2);      // slice starting at element 1, length 5, stride 2

   TRACE(EvenElements)(U[EvenElements]);

   U[EvenElements]*= 10;                      // doesn't work
   U[EvenElements] = U[EvenElements]*10;      // doesn't work
   U[EvenElements] = Vector<double>(5,10.);     // works
   //   U[EvenElements] = 10;                      // doesn't work

   U[EvenElements] = FixedVector<double>(size(EvenElements), 10);

   fill(U[EvenElements], 5);

   TRACE(U);                                    // here, U = (10,1,10,1,10,1,10,1,10,1)

}
