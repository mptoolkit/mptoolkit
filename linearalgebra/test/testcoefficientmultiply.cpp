// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testcoefficientmultiply.cpp
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

#include "linearalgebra/matrix.h"
#include "linearalgebra/coefficient_operations.h"
#include "linearalgebra/matrixtransform.h"
#include "linearalgebra/matrixtranspose.h"
#include "linearalgebra/matrixaddition.h"
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/hashvector.h"

using namespace LinearAlgebra;

struct CF
{
   double operator()(int i, int k, int j) const
   {
      return i+j+k;
   }

   typedef double result_type;
};

template <typename Mtype, typename Ntype, typename Ctype>
void CheckCF()
{
   Mtype M(3,3);
   zero_all(M);
   M(0,0) = 1;
   M(0,1) = 5;
   M(0,2) = 4;
   M(1,0) = 3;
   M(1,2) = 7;
   M(2,0) = 9;
   M(2,1) = 6;
   M(2,2) = 2;

   Ntype N(3,3);
   zero_all(N);
   N(0,0) = 1;
   N(0,1) = 5;
   N(0,2) = 4;
   N(1,0) = 3;
   N(1,2) = 7;
   N(2,0) = 9;
   N(2,1) = 6;
   N(2,2) = 2;

   CHECK_CLOSE(M,N);

   Matrix<double> MMCFCheck(3,3);
   MMCFCheck(0,0) = 87;
   MMCFCheck(0,1) = 77;
   MMCFCheck(0,2) = 145;
   MMCFCheck(1,0) = 192;
   MMCFCheck(1,1) = 198;
   MMCFCheck(1,2) = 106;
   MMCFCheck(2,0) = 144;
   MMCFCheck(2,1) = 195;
   MMCFCheck(2,2) = 378;

   Ctype MMCF = coefficient_multiply(M, N, CF());
   CHECK_CLOSE(MMCF, MMCFCheck);

   CHECK_CLOSE(MMCF, coefficient_multiply(M, N, CF()));

   MMCF += coefficient_multiply(M, N, CF());
   MMCF *= 0.5;
   CHECK_CLOSE(MMCF, MMCFCheck);

   MMCF -= coefficient_multiply(M, N, CF());
   CHECK_CLOSE(norm_frob(MMCF), 0);
}

int main()
{
   typedef SparseMatrix<double, RowMajor> t1;
   typedef SparseMatrix<double, ColMajor> t2;

   typedef SparseMatrix<double, RowMajor, HashVector<double> > t3;
   typedef SparseMatrix<double, ColMajor, HashVector<double> > t4;

   typedef Matrix<double, RowMajor> t5;
   typedef Matrix<double, ColMajor> t6;

   CheckCF<t1,t1,t1>();
   CheckCF<t2,t1,t1>();
   CheckCF<t3,t1,t1>();
   CheckCF<t4,t1,t1>();

   CheckCF<t1,t2,t1>();
   CheckCF<t2,t2,t1>();
   CheckCF<t3,t2,t1>();
   CheckCF<t4,t2,t1>();

   CheckCF<t1,t3,t1>();
   CheckCF<t2,t3,t1>();
   CheckCF<t3,t3,t1>();
   CheckCF<t4,t3,t1>();

   CheckCF<t1,t4,t1>();
   CheckCF<t2,t4,t1>();
   CheckCF<t3,t4,t1>();
   CheckCF<t4,t4,t1>();


   CheckCF<t1,t1,t2>();
   CheckCF<t2,t1,t2>();
   CheckCF<t3,t1,t2>();
   CheckCF<t4,t1,t2>();

   CheckCF<t1,t2,t2>();
   CheckCF<t2,t2,t2>();
   CheckCF<t3,t2,t2>();
   CheckCF<t4,t2,t2>();

   CheckCF<t1,t3,t2>();
   CheckCF<t2,t3,t2>();
   CheckCF<t3,t3,t2>();
   CheckCF<t4,t3,t2>();

   CheckCF<t1,t4,t2>();
   CheckCF<t2,t4,t2>();
   CheckCF<t3,t4,t2>();
   CheckCF<t4,t4,t2>();


   CheckCF<t1,t1,t3>();
   CheckCF<t2,t1,t3>();
   CheckCF<t3,t1,t3>();
   CheckCF<t4,t1,t3>();

   CheckCF<t1,t2,t3>();
   CheckCF<t2,t2,t3>();
   CheckCF<t3,t2,t3>();
   CheckCF<t4,t2,t3>();

   CheckCF<t1,t3,t3>();
   CheckCF<t2,t3,t3>();
   CheckCF<t3,t3,t3>();
   CheckCF<t4,t3,t3>();

   CheckCF<t1,t4,t3>();
   CheckCF<t2,t4,t3>();
   CheckCF<t3,t4,t3>();
   CheckCF<t4,t4,t3>();


   CheckCF<t1,t1,t4>();
   CheckCF<t2,t1,t4>();
   CheckCF<t3,t1,t4>();
   CheckCF<t4,t1,t4>();

   CheckCF<t1,t2,t4>();
   CheckCF<t2,t2,t4>();
   CheckCF<t3,t2,t4>();
   CheckCF<t4,t2,t4>();

   CheckCF<t1,t3,t4>();
   CheckCF<t2,t3,t4>();
   CheckCF<t3,t3,t4>();
   CheckCF<t4,t3,t4>();

   CheckCF<t1,t4,t4>();
   CheckCF<t2,t4,t4>();
   CheckCF<t3,t4,t4>();
   CheckCF<t4,t4,t4>();


   CheckCF<t1,t1,t5>();
   CheckCF<t2,t1,t5>();
   CheckCF<t3,t1,t5>();
   CheckCF<t4,t1,t5>();
   CheckCF<t5,t1,t5>();
   CheckCF<t6,t1,t5>();

   CheckCF<t1,t2,t5>();
   CheckCF<t2,t2,t5>();
   CheckCF<t3,t2,t5>();
   CheckCF<t4,t2,t5>();
   CheckCF<t5,t2,t5>();
   CheckCF<t6,t2,t5>();

   CheckCF<t1,t3,t5>();
   CheckCF<t2,t3,t5>();
   CheckCF<t3,t3,t5>();
   CheckCF<t4,t3,t5>();
   CheckCF<t5,t3,t5>();
   CheckCF<t6,t3,t5>();

   CheckCF<t1,t4,t5>();
   CheckCF<t2,t4,t5>();
   CheckCF<t3,t4,t5>();
   CheckCF<t4,t4,t5>();
   CheckCF<t5,t4,t5>();
   CheckCF<t6,t4,t5>();

   CheckCF<t1,t5,t5>();
   CheckCF<t2,t5,t5>();
   CheckCF<t3,t5,t5>();
   CheckCF<t4,t5,t5>();
   CheckCF<t5,t5,t5>();
   CheckCF<t6,t5,t5>();

   CheckCF<t1,t6,t5>();
   CheckCF<t2,t6,t5>();
   CheckCF<t3,t6,t5>();
   CheckCF<t4,t6,t5>();
   CheckCF<t5,t6,t5>();
   CheckCF<t6,t6,t5>();


   CheckCF<t1,t1,t6>();
   CheckCF<t2,t1,t6>();
   CheckCF<t3,t1,t6>();
   CheckCF<t4,t1,t6>();
   CheckCF<t5,t1,t6>();
   CheckCF<t6,t1,t6>();

   CheckCF<t1,t2,t6>();
   CheckCF<t2,t2,t6>();
   CheckCF<t3,t2,t6>();
   CheckCF<t4,t2,t6>();
   CheckCF<t5,t2,t6>();
   CheckCF<t6,t2,t6>();

   CheckCF<t1,t3,t6>();
   CheckCF<t2,t3,t6>();
   CheckCF<t3,t3,t6>();
   CheckCF<t4,t3,t6>();
   CheckCF<t5,t3,t6>();
   CheckCF<t6,t3,t6>();

   CheckCF<t1,t4,t6>();
   CheckCF<t2,t4,t6>();
   CheckCF<t3,t4,t6>();
   CheckCF<t4,t4,t6>();
   CheckCF<t5,t4,t6>();
   CheckCF<t6,t4,t6>();

   CheckCF<t1,t5,t6>();
   CheckCF<t2,t5,t6>();
   CheckCF<t3,t5,t6>();
   CheckCF<t4,t5,t6>();
   CheckCF<t5,t5,t6>();
   CheckCF<t6,t5,t6>();

   CheckCF<t1,t6,t6>();
   CheckCF<t2,t6,t6>();
   CheckCF<t3,t6,t6>();
   CheckCF<t4,t6,t6>();
   CheckCF<t5,t6,t6>();
   CheckCF<t6,t6,t6>();
}
