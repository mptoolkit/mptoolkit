// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testlinearsolvehpd.cpp
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

#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/vector.h"
#include "linearalgebra/matrixtranspose.h"
#include "linearalgebra/matrixmatrixmultiplication.h"
#include "linearalgebra/matrixaddition.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/matrixtransform.h"
#include <complex>

using namespace LinearAlgebra;

int main()
{
   Matrix<std::complex<double> > M(10,10);  // left hand side
   Matrix<std::complex<double> > R(10,2);   // 2 right hand side vectors

   for (int i = 0; i < 10; ++i)
   {
      for (int j = 0; j < 10; ++j)
      {
         M(i,j) = std::complex<double>(i * j + i + j, i-j);
         if (i == j) M(i,j) += 1000.0;
      }

      for (int j = 0; j < 2; ++j)
      {
         R(i,j) = std::complex<double>(j+i, j+2 * i);
      }
   }

   Matrix<std::complex<double> > S = LinearSolveHPD(M,R);
   Matrix<std::complex<double> > RCheck = M * S;

   CHECK_CLOSE(R, RCheck);

   CHECK_CLOSE(real(R), M * LinearSolveHPD(M, real(R)));
   CHECK_CLOSE(R, real(M) * LinearSolveHPD(real(M), R));
   CHECK_CLOSE(real(R), real(M) * LinearSolveHPD(real(M), real(R)));
}
