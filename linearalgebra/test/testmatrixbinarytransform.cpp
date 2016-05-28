// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testmatrixbinarytransform.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   LinearAlgebra::Matrix<double> M(3,3, 2.0);
   LinearAlgebra::Matrix<double> N(3,3, 3.0);
   
   LinearAlgebra::Matrix<double> R;

   R = transform(M, N, LinearAlgebra::Multiplication<double, double>());

   TRACE(R);

   LinearAlgebra::Matrix<double, RowMajor> MM(2,3, 2.0);
   LinearAlgebra::Matrix<double, ColMajor> NN(2,3, 3.0);

   TRACE(transform(MM, NN, LinearAlgebra::Multiplication<double, double>()));

   TRACE(MM+NN);

   TRACE(element_prod(MM,NN));
}
