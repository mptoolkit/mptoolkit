// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testvector.cpp
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
#include "testdensevector.h"

int main()
{
   test_dense_vector<LinearAlgebra::Vector<double> >();

   test_real_vector<LinearAlgebra::Vector<double> >();

   typedef LinearAlgebra::Vector<LinearAlgebra::Vector<double> > vt;

   test_ctor<vt>();

   LinearAlgebra::Vector<double> v1(4,1);
   test_dense_single<vt>(v1);

   vt v2(10, v1);
   test_assign(v2);

   test_double_negate<vt>();
   test_real_scalar_vector_nop<vt>(v2);
   test_equal_to(v2);

   LinearAlgebra::Vector<double> V(1);
   V = V[LinearAlgebra::Range(0,0)];
   CHECK_EQUAL(size(V), 0);
}
