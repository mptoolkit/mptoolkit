// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testvectorrange.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/vector.h"

using namespace LinearAlgebra;

int main()
{
   Vector<int> v = range(0,5);
   TRACE(v);
   v = v[range(2,5)];
   TRACE(v);
}
