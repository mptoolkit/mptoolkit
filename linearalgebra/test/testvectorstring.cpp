// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testvectorstring.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/vector.h"
#include <string>

using namespace LinearAlgebra;

int main()
{
   Vector<std::string> X(3);
   X[0] = "hello";
   X[1] = "world";
   X[2] = "goodbye";

   X = X[range(1,2)];
   CHECK_EQUAL(size(X), 1);
   CHECK_EQUAL(X[0], "world");
}
