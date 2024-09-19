// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testmatrixserialize.cpp
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
#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/pstreamio.h"
#include "pstream/pfilestream.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double> M = random_matrix<double>(10,20);

   {
      PStream::opfilestream Out;
      Out.open("hello.test.xdr", O_WRONLY | O_TRUNC | O_CREAT);
      Out.set_format(PStream::format::XDR);
      Out << M;
      Out.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Matrix<double> MCheck;
      In >> MCheck;
      CHECK_CLOSE(M, MCheck);
      In.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Matrix<double, ColMajor> MCheck;
      In >> MCheck;
      CHECK_CLOSE(M, MCheck);
      In.close();
   }

   SparseMatrix<double> N(10,20);
   N(0,1) = 5;
   N(1,2) = -10;

   {
      PStream::opfilestream Out;
      Out.open("hello.test.xdr", O_WRONLY | O_TRUNC | O_CREAT);
      Out.set_format(PStream::format::XDR);
      Out << N;
      Out.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      SparseMatrix<double, RowMajor> NCheck;
      In >> NCheck;
      CHECK_CLOSE(N, NCheck);
      In.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      SparseMatrix<double, ColMajor> NCheck;
      In >> NCheck;
      CHECK_CLOSE(N, NCheck);
      In.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Matrix<double, RowMajor> NCheck;
      In >> NCheck;
      CHECK_CLOSE(N, NCheck);
      In.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Matrix<double, ColMajor> NCheck;
      In >> NCheck;
      CHECK_CLOSE(N, NCheck);
      In.close();
   }
}
