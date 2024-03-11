// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/testdgemm.cpp
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

#include "densematrix.h"

int main()
{
   DenseMatrix<double> Result(4,5,0.0);
   DenseMatrix<double> M1(4,6,1.0);
   DenseMatrix<double> M2(6,5,2.0);
   for (int i = 0; i < M1.rows(); ++i)
   {
      for (int j = 0; j < M1.cols(); ++j)
      {
         M1(i,j) = i-j;
      }
   }

   for (int i = 0; i < M2.rows(); ++i)
   {
      for (int j = 0; j < M2.cols(); ++j)
      {
         M2(i,j) = i*j;
      }
   }


   add_product_AB(Result, M1, M2);

   BLAS::DGS::DebugPrintDgemm();
}
