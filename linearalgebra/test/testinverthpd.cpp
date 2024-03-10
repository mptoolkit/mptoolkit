// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testinverthpd.cpp
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

int main()
{
   int const dim = 8;

   Matrix<std::complex<double> > I = diagonal_matrix(FixedVector<double>(dim, 1.0));
   Matrix<std::complex<double> > H = random_matrix<std::complex<double> >(dim, dim);
   H = H * herm(H) + I;

   Matrix<std::complex<double> > HInv(H);

   InvertHPD(HInv);

   CHECK_CLOSE(HInv * H, I)(H)(HInv);
   CHECK_CLOSE(H * HInv, I)(H)(HInv);
}
