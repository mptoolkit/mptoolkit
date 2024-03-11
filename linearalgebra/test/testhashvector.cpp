// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testhashvector.cpp
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

#include "linearalgebra/hashvector.h"
#include "testvectorsparse.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrixmatrixmultiplication.h"

int main()
{
   test_real_sparse_vector<LinearAlgebra::HashVector<double> >();

#if 0
   LinearAlgebra::HashVector<LinearAlgebra::Matrix<double> > v1(10), v2(10);
   TRACE(typeid(parallel_prod(v1,v2)).name());
#endif
}
