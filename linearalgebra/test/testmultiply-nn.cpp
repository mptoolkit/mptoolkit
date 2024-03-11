// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testmultiply-nn.cpp
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

#include "linearalgebra/matrix.h"
#include "linearalgebra/matrixmatrixmultiplication.h"
#include "linearalgebra/matrixtransform.h"
#include "linearalgebra/matrixtranspose.h"
#include "linearalgebra/matrixaddition.h"
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/hashvector.h"

using namespace LinearAlgebra;

template <typename Mtype, typename Ntype, typename Ctype>
void CheckMult1()
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

   //   TRACE(norm_frob(M))(norm_frob(N));

   //   TRACE(norm_frob(M-N));

   CHECK_CLOSE(M,N);

   Matrix<double> MMCFCheck(3,3);
   MMCFCheck(0,0) = 52;
   MMCFCheck(0,1) = 29;
   MMCFCheck(0,2) = 47;
   MMCFCheck(1,0) = 66;
   MMCFCheck(1,1) = 57;
   MMCFCheck(1,2) = 26;
   MMCFCheck(2,0) = 45;
   MMCFCheck(2,1) = 57;
   MMCFCheck(2,2) = 82;

   Ctype MMCF = M * N;

   CHECK_CLOSE(MMCF, MMCFCheck);
   CHECK_CLOSE(MMCF, M*N);

   MMCF += M*N;
   MMCF *= 0.5;
   CHECK_CLOSE(MMCF, MMCFCheck);

   MMCF -= 1.0 * M*N;
   CHECK_CLOSE(norm_frob(MMCF), 0);
}

template <typename Mtype, typename Ntype, typename Ctype>
void CheckMult2()
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
   MMCFCheck(0,0) = 52;
   MMCFCheck(0,1) = 29;
   MMCFCheck(0,2) = 47;
   MMCFCheck(1,0) = 66;
   MMCFCheck(1,1) = 57;
   MMCFCheck(1,2) = 26;
   MMCFCheck(2,0) = 45;
   MMCFCheck(2,1) = 57;
   MMCFCheck(2,2) = 82;

   Ctype MMCF = M * N;

   CHECK_CLOSE(MMCF, MMCFCheck);
   CHECK_CLOSE(MMCF, M*N);

   MMCF += M*N;
   MMCF *= 0.5;
   CHECK_CLOSE(MMCF, MMCFCheck);

   MMCF -= 1.0 * M*N;
   CHECK_CLOSE(norm_frob(MMCF), 0);
}

template <typename Mtype, typename Ntype, typename Ctype>
void CheckMult()
{
   CheckMult1<Mtype, Ntype, Ctype>();
   CheckMult2<Mtype, Ntype, Ctype>();
   //   CheckMult3<Mtype, Ntype, Ctype>();
}

int main()
{
   typedef SparseMatrix<double, RowMajor> t1;
   typedef SparseMatrix<double, ColMajor> t2;

#if 0
   typedef t1 t3;
   typedef t2 t4;
#else
   typedef SparseMatrix<double, RowMajor, HashVector<double> > t3;
   typedef SparseMatrix<double, ColMajor, HashVector<double> > t4;
#endif

   typedef Matrix<double, RowMajor> t5;
   typedef Matrix<double, ColMajor> t6;

   CheckMult<t1,t1,t1>();
   CheckMult<t2,t1,t1>();
   CheckMult<t3,t1,t1>();
   CheckMult<t4,t1,t1>();

   CheckMult<t1,t2,t1>();
   CheckMult<t2,t2,t1>();
   CheckMult<t3,t2,t1>();
   CheckMult<t4,t2,t1>();

   CheckMult<t1,t3,t1>();
   CheckMult<t2,t3,t1>();
   CheckMult<t3,t3,t1>();
   CheckMult<t4,t3,t1>();

   CheckMult<t1,t4,t1>();
   CheckMult<t2,t4,t1>();
   CheckMult<t3,t4,t1>();
   CheckMult<t4,t4,t1>();


   CheckMult<t1,t1,t2>();
   CheckMult<t2,t1,t2>();
   CheckMult<t3,t1,t2>();
   CheckMult<t4,t1,t2>();

   CheckMult<t1,t2,t2>();
   CheckMult<t2,t2,t2>();
   CheckMult<t3,t2,t2>();
   CheckMult<t4,t2,t2>();

   CheckMult<t1,t3,t2>();
   CheckMult<t2,t3,t2>();
   CheckMult<t3,t3,t2>();
   CheckMult<t4,t3,t2>();

   CheckMult<t1,t4,t2>();
   CheckMult<t2,t4,t2>();
   CheckMult<t3,t4,t2>();
   CheckMult<t4,t4,t2>();


   CheckMult<t1,t1,t3>();
   CheckMult<t2,t1,t3>();
   CheckMult<t3,t1,t3>();
   CheckMult<t4,t1,t3>();

   CheckMult<t1,t2,t3>();
   CheckMult<t2,t2,t3>();
   CheckMult<t3,t2,t3>();
   CheckMult<t4,t2,t3>();

   CheckMult<t1,t3,t3>();
   CheckMult<t2,t3,t3>();
   CheckMult<t3,t3,t3>();
   CheckMult<t4,t3,t3>();

   CheckMult<t1,t4,t3>();
   CheckMult<t2,t4,t3>();
   CheckMult<t3,t4,t3>();
   CheckMult<t4,t4,t3>();


   CheckMult<t1,t1,t4>();
   CheckMult<t2,t1,t4>();
   CheckMult<t3,t1,t4>();
   CheckMult<t4,t1,t4>();

   CheckMult<t1,t2,t4>();
   CheckMult<t2,t2,t4>();
   CheckMult<t3,t2,t4>();
   CheckMult<t4,t2,t4>();

   CheckMult<t1,t3,t4>();
   CheckMult<t2,t3,t4>();
   CheckMult<t3,t3,t4>();
   CheckMult<t4,t3,t4>();

   CheckMult<t1,t4,t4>();
   CheckMult<t2,t4,t4>();
   CheckMult<t3,t4,t4>();
   CheckMult<t4,t4,t4>();


   CheckMult<t1,t1,t5>();
   CheckMult<t2,t1,t5>();
   CheckMult<t3,t1,t5>();
   CheckMult<t4,t1,t5>();
   CheckMult<t5,t1,t5>();
   CheckMult<t6,t1,t5>();

   CheckMult<t1,t2,t5>();
   CheckMult<t2,t2,t5>();
   CheckMult<t3,t2,t5>();
   CheckMult<t4,t2,t5>();
   CheckMult<t5,t2,t5>();
   CheckMult<t6,t2,t5>();

   CheckMult<t1,t3,t5>();
   CheckMult<t2,t3,t5>();
   CheckMult<t3,t3,t5>();
   CheckMult<t4,t3,t5>();
   CheckMult<t5,t3,t5>();
   CheckMult<t6,t3,t5>();

   CheckMult<t1,t4,t5>();
   CheckMult<t2,t4,t5>();
   CheckMult<t3,t4,t5>();
   CheckMult<t4,t4,t5>();
   CheckMult<t5,t4,t5>();
   CheckMult<t6,t4,t5>();

   CheckMult<t1,t5,t5>();
   CheckMult<t2,t5,t5>();
   CheckMult<t3,t5,t5>();
   CheckMult<t4,t5,t5>();
   CheckMult<t5,t5,t5>();
   CheckMult<t6,t5,t5>();

   CheckMult<t1,t6,t5>();
   CheckMult<t2,t6,t5>();
   CheckMult<t3,t6,t5>();
   CheckMult<t4,t6,t5>();
   CheckMult<t5,t6,t5>();
   CheckMult<t6,t6,t5>();


   CheckMult<t1,t1,t6>();
   CheckMult<t2,t1,t6>();
   CheckMult<t3,t1,t6>();
   CheckMult<t4,t1,t6>();
   CheckMult<t5,t1,t6>();
   CheckMult<t6,t1,t6>();

   CheckMult<t1,t2,t6>();
   CheckMult<t2,t2,t6>();
   CheckMult<t3,t2,t6>();
   CheckMult<t4,t2,t6>();
   CheckMult<t5,t2,t6>();
   CheckMult<t6,t2,t6>();

   CheckMult<t1,t3,t6>();
   CheckMult<t2,t3,t6>();
   CheckMult<t3,t3,t6>();
   CheckMult<t4,t3,t6>();
   CheckMult<t5,t3,t6>();
   CheckMult<t6,t3,t6>();

   CheckMult<t1,t4,t6>();
   CheckMult<t2,t4,t6>();
   CheckMult<t3,t4,t6>();
   CheckMult<t4,t4,t6>();
   CheckMult<t5,t4,t6>();
   CheckMult<t6,t4,t6>();

   CheckMult<t1,t5,t6>();
   CheckMult<t2,t5,t6>();
   CheckMult<t3,t5,t6>();
   CheckMult<t4,t5,t6>();
   CheckMult<t5,t5,t6>();
   CheckMult<t6,t5,t6>();

   CheckMult<t1,t6,t6>();
   CheckMult<t2,t6,t6>();
   CheckMult<t3,t6,t6>();
   CheckMult<t4,t6,t6>();
   CheckMult<t5,t6,t6>();
   CheckMult<t6,t6,t6>();
}
