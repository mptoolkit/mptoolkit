// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/testdgemm.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
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
