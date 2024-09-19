// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testinverttri.cpp
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
#include "linearalgebra/fixedvector.h"

#include <iostream>
#include <iomanip>

using namespace LinearAlgebra;
typedef std::complex<double> complex;

template <typename Mat>
void ZeroUpper(Mat& M)
{
   for (unsigned i = 0; i < M.size1(); ++i)
   {
      for (unsigned j = i+1; j < M.size2(); ++j)
      {
         M(i,j) = 0.0;
      }
   }
}

template <typename Mat>
void ZeroLower(Mat& M)
{
   for (unsigned j = 0; j < M.size2(); ++j)
   {
      for (unsigned i = j+1; i < M.size1(); ++i)
      {
         M(i,j) = 0.0;
      }
   }
}

int main()
{
   int const dim = 8;

   Matrix<std::complex<double> > I = diagonal_matrix(FixedVector<double>(dim, 1.0));
   Matrix<std::complex<double> > H = random_matrix<std::complex<double> >(dim, dim);
   H = H * herm(H) + I;
   CholeskyFactorizeUpper(H);
   ZeroLower(H);

   Matrix<std::complex<double> > HInv(H);
   InvertUpperTriangular(HInv);
   CHECK_CLOSE(HInv * H, I)(H)(HInv);
   CHECK_CLOSE(H * HInv, I)(H)(HInv);

   Matrix<std::complex<double> > H2 = herm(H);
   Matrix<std::complex<double> > Hi2 = H2;
   InvertLowerTriangular(Hi2);
   CHECK_CLOSE(Hi2 * H2, I)(H2)(Hi2);
   CHECK_CLOSE(H2 * Hi2, I)(H2)(Hi2);

   Matrix<complex, ColMajor> H3 = H;
   Matrix<std::complex<double> > Hi3 = H3;
   InvertUpperTriangular(Hi3);
   CHECK_CLOSE(Hi3 * H3, I)(H3)(Hi3);
   CHECK_CLOSE(H3 * Hi3, I)(H3)(Hi3);

   Matrix<complex, ColMajor> H4 = herm(H);
   Matrix<std::complex<double> > Hi4 = H4;
   InvertLowerTriangular(Hi4);
   CHECK_CLOSE(Hi4 * H4, I)(H4)(Hi4);
   CHECK_CLOSE(H4 * Hi4, I)(H4)(Hi4);
}
