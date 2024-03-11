// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testvectorstring.cpp
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
