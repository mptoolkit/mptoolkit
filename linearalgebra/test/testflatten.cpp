// -*- C++ -*- $Id$

#include "linearalgebra/matrix.h"

using namespace LinearAlgebra;

int main()
{
   Matrix<double, RowMajor> M(3,4,0);
   M(0,0) = 1;
   M(0,1) = 2;
   M(0,2) = 3;
   M(0,3) = 4;
   M(1,0) = 11;
   M(1,1) = 12;
   M(1,2) = 13;
   M(1,3) = 14;
   M(2,0) = 21;
   M(2,1) = 22;
   M(2,2) = 23;
   M(2,3) = 24;

   Vector<double> V(12,0);
   V[0] = 1;
   V[1] = 2;
   V[2] = 3;
   V[3] = 4;
   V[4] = 11;
   V[5] = 12;
   V[6] = 13;
   V[7] = 14;
   V[8] = 21;
   V[9] = 22;
   V[10] = 23;
   V[11] = 24;

   CHECK_EQUAL(flatten_rows(M), V);
   CHECK_EQUAL(flatten_cols(transpose(M)), V);
}
