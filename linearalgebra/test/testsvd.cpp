// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// linearalgebra/test/testsvd.cpp
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

#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix.h"
#include "linearalgebra/diagonalmatrix.h"
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

int main()
{
   for (int i = 1; i <= 10; ++i)
   {
      for (int j = 1; j <= 10; ++j)
      {
	 Matrix<std::complex<double>> M = random_matrix<std::complex<double>>(i,j);

	 Matrix<std::complex<double>> U, Vh;
	 DiagonalMatrix<double> D;
	 
	 SingularValueDecompositionFull(M, U, D, Vh);
	 
	 CHECK(norm_frob(U*D*Vh - M) < 1E-10)(M)(U)(D)(Vh);
      }
   }
}
