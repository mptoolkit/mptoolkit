// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testsum.cpp
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
#include "linearalgebra/matrix.h"
#include "linearalgebra/scalar.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> v = Range(1,10);
   CHECK_EQUAL(sum(v), 45.0);

   Matrix<double> M(3,3,1.0);
   Vector<Matrix<double> > vm(3, M);

   Matrix<double> N(3,3,0.0);

   typedef boost::mpl::print<Sum<Vector<Matrix<double> > >::result_type>::type d;


   N += sum(vm);
   CHECK_EQUAL(N, 3*M);

   Vector<Matrix<double> > vn;
   N += sum(vn);
   CHECK_EQUAL(N, 3*M);
}
