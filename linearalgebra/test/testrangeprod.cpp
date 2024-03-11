// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// linearalgebra/test/testrangeprod.cpp
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
#include "linearalgebra/matrix_utility.h"

using namespace LinearAlgebra;

using tracer::typeid_name;

int main()
{
   Matrix<double> M(random_matrix<double>(40, 40));

   Matrix<double> R1(M(range(10,20), all));
   Matrix<double> R1T(transpose(R1));
   Matrix<double> R1R1TCheck = R1 * R1T;
   Matrix<double> R1R1T = M(range(10,20), all) * transpose(M(range(10,20), all));
   CHECK_CLOSE(R1R1T, R1R1TCheck);
   R1R1T = M(range(10,20), all) * conj(transpose(M(range(10,20), all)));
   CHECK_CLOSE(R1R1T, R1R1TCheck);
   R1R1T = conj(M(range(10,20), all)) * conj(transpose(M(range(10,20), all)));
   CHECK_CLOSE(R1R1T, R1R1TCheck);

   CHECK_CLOSE(transpose(M(all, range(10,20))), transpose(M)(range(10, 20), all))(M);

   Matrix<double> M1(M(range(10,20), range(20,30)));
   Matrix<double> M2(M(range(20,30), range(10,20)));

   Matrix<double> M1M2Check = M1 * M2;
   Matrix<double> M1M2 = M(range(10,20), range(20,30)) * M(range(20,30), range(10,20));
   CHECK_CLOSE(M1M2, M1M2Check);

   Matrix<double> A = random_matrix<double>(10,10);
   TRACE(A);
   A(all, range(1,2)) = A(all, range(1,2)) * random_matrix<double>(1, 1);
   TRACE(A);
}
