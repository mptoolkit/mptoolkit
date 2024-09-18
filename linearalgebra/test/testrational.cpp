// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testrational.cpp
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

#include "linearalgebra/vector.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/gmpalgebra.h"

using namespace LinearAlgebra;

int main()
{
   Vector<gmp::rational> Vec1(10, gmp::rational(0.1));
   TRACE(Vec1);

   Vector<gmp::rational> Vec2(10, gmp::rational(2.0));

   Vector<gmp::rational> Vec3 = Vec1 + Vec2;

   TRACE(norm_frob_sq(Vec3));

   TRACE(Vec3);

   TRACE(inner_prod(Vec1, Vec2));
}
