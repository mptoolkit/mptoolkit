// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testqr.cpp
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

#include "linearalgebra/matrix.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/fixedvector.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;

typedef std::complex<double> complex;

void Test(int m, int n)
{
   {
      Matrix<std::complex<double>> H = random_matrix<std::complex<double>>(m, n);
      Matrix<std::complex<double>> R(H);
      //TRACE(H);
      Matrix<std::complex<double>> Q = QR_FactorizeFull(R);
      CHECK_CLOSE(H,Matrix<complex>(Q*R))(Q)(R)(conj(Q)*R)(herm(Q)*R)(trans(Q)*R);
   }

   {
      Matrix<double> H = random_matrix<double>(m, n);
      Matrix<double> R(H);
      //TRACE(H);
      Matrix<double> Q = QR_FactorizeFull(R);
      CHECK_CLOSE(H,Matrix<double>(Q*R))(Q)(R)(conj(Q)*R)(herm(Q)*R)(trans(Q)*R);
   }
}

int main()
{
   Test(2,2);
   Test(3,3);
   Test(4,4);
   Test(2,1);
   Test(2,3);
   Test(10,5);
   Test(5,10);
}
