// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testscalarmatrix.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "linearalgebra/scalarmatrix.h"
#include "linearalgebra/scalar.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/matrixsection.h"
#include "linearalgebra/sparsematrix.h"

#include "linearalgebra/matrixbinarytransform.h"

using namespace LinearAlgebra;

using tracer::typeid_name;

int main()
{
   ScalarMatrix<double> M(3,3,2.0);

   ScalarMatrix<double> N(3,3,5.0);

   TRACE(M);

   TRACE(M*N);
   TRACE(M+N);

   TRACE(typeid_name(M+N));

   CHECK_CLOSE(M*N, N+N);

   CHECK_CLOSE(M*N, N*2.0);
   CHECK_CLOSE(M*N, 2.0*N);
   CHECK_CLOSE(M*N, N*2);
   CHECK_CLOSE(M*N, 2*N);

   for (size_type i = 0; i < size1(M); ++i)
   {
      for (size_type j = 0; j < size2(M); ++j)
      {
         if (i == j)
            CHECK_EQUAL(M(i,j), 2);
         else
            CHECK_EQUAL(M(i,j), 0);
      }
   }

   Matrix<double> Test(3,3,0);
   Test(0,0) = 5;
   Test(0,1) = 2;
   Test(1,0) = 4;
   Test(1,1) = 8;
   Test(2,1) = 9;
   Test(2,2) = 10;

   TRACE(M+Test);

   TRACE(Test+M);

   TRACE(M * Test);

   TRACE(typeid_name(M*Test));
   TRACE(typeid_name(Test*M));

   TRACE(typeid_name(M+Test));
   TRACE(typeid_name(Test+M));

   TRACE(typeid_name(Test-M));

   TRACE(Test-M);
   TRACE(M-Test);

   Test += M;

   ScalarMatrix<std::complex<double> > C(3,3,std::complex<double>(1.0,2.0));
   real(C) *= 4.0;
   TRACE(C);

   TRACE(direct_product(M, Test));
}
