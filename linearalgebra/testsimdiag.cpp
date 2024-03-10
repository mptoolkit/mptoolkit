// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/testsimdiag.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
// Copyright (C) 2022 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "simdiag.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<std::complex<double>> A(3,3, 0.0);
   Matrix<std::complex<double>> B(3,3, 0.0);

   A(0,0) = 1.0;
   A(1,1) = 1.0;
   A(2,2) = 1.0;
   A(1,2) = {0.0,0.5};
   A(2,1) = {0.0,-0.5};

   B(0,0) = 2.0;
   B(1,1) = 4.0;
   B(2,2) = 4.0;

   std::vector<Matrix<std::complex<double>>> M = {A,B};

   LinearAlgebra::Matrix<std::complex<double>> Q;
   std::vector<LinearAlgebra::Vector<double>> EVal;

   std::tie(Q, EVal) = SimultaneousDiagonalizeHermitian(M);

   TRACE(Q);

   TRACE(EVal[0]);
   TRACE(herm(Q)*A*Q);

   TRACE(EVal[1]);
   TRACE(herm(Q)*B*Q);
}
