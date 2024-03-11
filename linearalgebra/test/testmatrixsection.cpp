// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testmatrixsection.cpp
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
#include "linearalgebra/matrixsection.h"
#include "linearalgebra/matrixtranspose.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M(3,3,0);
   M(0,0) = 1;
   M(0,1) = 5;
   M(0,2) = 4;
   M(1,0) = 3;
   M(1,1) = 8;
   M(1,2) = 7;
   M(2,0) = 9;
   M(2,1) = 6;
   M(2,2) = 2;

   TRACE(M);

   Range r(0,2);
   Slice s(0,2,2);


   BOOST_MPL_ASSERT((is_proxy_reference<Range const&>));
   BOOST_MPL_ASSERT((is_const_proxy_reference<Range const&>));

   MatrixSection<Matrix<double>&, Range const&, Slice const&> MySec(M, r, s);
   TRACE(MySec);

   get_element(MySec, 1,1) = -10;
   TRACE(MySec);
   TRACE(M);

   Matrix<double> Test(MySec);

   MySec = transpose(Test);

   TRACE(Test)(MySec)(M);

   TRACE(project(M, r, s));

   TRACE(M);
   M(r,s) *= 2.0;
   TRACE(M);
}
