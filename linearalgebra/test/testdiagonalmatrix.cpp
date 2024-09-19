// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testdiagonalmatrix.cpp
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

#include "linearalgebra/diagonalmatrix.h"
#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   DiagonalMatrix<double> M(5,5);
   M.diagonal() = range(1,6);

   TRACE(M);

   Matrix<double> Mat = M;
   Mat(0,1) = 1;
   Mat(1,0) = 2;
   Mat(0,2) = 3;
   Mat(0,4) = 7;
   TRACE(Mat);

   Matrix<double> Mat2 = Mat*M;

   TRACE(Mat);
   //   TRACE(Matrix<double>(M*Mat));
   TRACE(Matrix<double>(Mat*M));

   TRACE(M*M);

   TRACE(norm_frob_sq(M*M));

   //   CHECK_EQUAL(norm_frob_sq(M*M), norm_frob_sq(Mat2));
}
