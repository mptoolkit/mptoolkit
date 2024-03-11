// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testdiagonalizehermitian.cpp
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
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

int main()
{
   int const dim = 8;
   Matrix<std::complex<double> > H = random_matrix<std::complex<double> >(dim, dim);
   H=H+herm(H);

   // simplest case
   Matrix<std::complex<double> > HCopy(H);
   Vector<double> Eigen = DiagonalizeHermitian(H);
   CHECK_CLOSE(H * HCopy * herm(H), diagonal_matrix(Eigen));

   Matrix<std::complex<double> > T = H;
   H = HCopy;

   // column major matrix
   Matrix<std::complex<double>, ColMajor> HCol = H;
   Vector<double> Eigen2 = DiagonalizeHermitian(HCol);
   CHECK_CLOSE(Eigen2, Eigen);
   CHECK_CLOSE(HCol, T);
   // a slice expression
   Vector<double> EigenSlice = DiagonalizeHermitian(H(Slice(2,3,2), Slice(2,3,2)));
   Matrix<std::complex<double> > HSlice = H(Slice(2,3,2), Slice(2,3,2));
   Matrix<std::complex<double> > HCopySlice = HCopy(Slice(2,3,2), Slice(2,3,2));
   CHECK_CLOSE(HSlice * HCopySlice * herm(HSlice), diagonal_matrix(EigenSlice));
}
